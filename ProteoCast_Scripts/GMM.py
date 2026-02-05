import numpy as np
import os
import seaborn as sns
import pandas as pd
from sklearn import mixture
import matplotlib.pyplot as plt
import pickle
import sys
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

class GMM():
    def __init__(self,
            GEMMEmat_location="",
            prot_name="",
            nb_clusters=2,
            coVar_type='full', 
            outPickle='',
            #nonconfres=[],
            plot=True,
            alphabet='ACDEFGHIKLMNPQRSTVWY'
            ):
            self.GEMMEmat_location = GEMMEmat_location
            self.alphabet = list(alphabet)
            self.prot_name=prot_name
            self.outputDirPickle = outPickle
            #self.NonConf_resi=nonconfres
            self.nb_clusters=nb_clusters
            self.coVar_type=coVar_type
            self.plot=plot
            #self.SingleMutationsVector()
            self.GMM()


    def GMM(self):

        if not os.path.exists(self.outputDirPickle):
            print(f"Warning: The directory GMM does not exist. Creating it.")
            os.makedirs(self.outputDirPickle)

        self.df_single_mutations = pd.read_csv(self.GEMMEmat_location)
        ## if the Global Confidence is high then we take scores at high local confidence; otherwise everything with a warning message
        if True in self.df_single_mutations['LocalConfidence'].unique():
            self.df_single_mutationsConf = self.df_single_mutations.loc[self.df_single_mutations['LocalConfidence']==True]
        else:
            self.df_single_mutationsConf = self.df_single_mutations
        self.df_single_mutationsConf = self.df_single_mutationsConf.loc[self.df_single_mutationsConf['Mutation'].str[0] != self.df_single_mutationsConf['Mutation'].str[-1]]
        self.df_single_mutationsConf = self.df_single_mutationsConf.reset_index(drop=True)
        self.X_train = np.array(self.df_single_mutationsConf['Variant_score']).reshape(-1, 1)
        self.model_GMM = mixture.GaussianMixture(n_components=self.nb_clusters, covariance_type=self.coVar_type,max_iter=1000, n_init=30,tol=1e-5)
        self.model_GMM.fit(self.X_train)

        labels=self.model_GMM.predict(self.X_train)
        temp_labels = np.copy(labels) 

        self.cluster_means = self.model_GMM.means_.flatten()
        min_mean_index = np.argmin(self.cluster_means)
        max_mean_index = np.argmax(self.cluster_means)
        
        temp_labels[labels == min_mean_index] = 2
        temp_labels[labels == max_mean_index] = 0
        temp_labels[(labels != min_mean_index) & (labels != max_mean_index)] = 1

        self.df_single_mutationsConf.loc[:, 'labels'] = temp_labels
                                            
        ####### THRESHOLDS ########

        component_share = self.model_GMM.predict_proba(self.X_train)  # Shape: (#values, #clusters)
        self.sorted_indices = np.argsort(self.cluster_means)
        self.cluster_means = self.cluster_means[self.sorted_indices]

        self.component_shareOrdered = component_share[:, self.sorted_indices]
        df = pd.DataFrame(self.component_shareOrdered, columns=['P1', 'P2', 'P3'])
        df = df.reset_index(drop=True)
        df = pd.concat([self.df_single_mutationsConf, df], axis=1)
        df_round = df.round(3)    
        df_round = df_round.loc[(df_round['Variant_score']>=self.cluster_means[0])&(df_round['Variant_score']<=self.cluster_means[2])]
        df_round.sort_values('Variant_score', ignore_index=True, inplace=True)
        proba_thresh=0.001

        try: 
            self.thresh = df_round.loc[(df_round['P1']>=proba_thresh)&(df_round['P2']>=proba_thresh)&(df_round['P1']<df_round['P2']), 'Variant_score'].tolist()[0]
        except:
            try:
                fr1 = df_round.loc[(df_round['P1']>=proba_thresh), 'Variant_score'].tolist()[-1]
            except:
                fr1 = df_round['Variant_score'].tolist()[0]
            fr2 = df_round.loc[(df_round['P2']>=proba_thresh), 'Variant_score'].tolist()[0]
            self.thresh  = (fr1+fr2)/2
        
        df_thresholds = pd.DataFrame([round(self.cluster_means[1],3), round(self.thresh,3)],  index = ['GMM_mild','GMM_impactful'],columns = ['Values'])

        df_thresholds.to_csv(f'{self.outputDirPickle}/{self.prot_name}_thresholds.csv')
        self.df_single_mutationsConf.to_csv(f'{self.outputDirPickle}/{self.prot_name}_labels.csv', index=False)
        filename = f'{self.outputDirPickle}/{self.prot_name}_GMM_model.pickle'
        pickle.dump(self.model_GMM, open(filename, 'wb'))
        if self.plot: 
            self.plotGMM()

        ## add Variant and Residue classification to ProteoCast.csv
        mut_scores = np.array(self.df_single_mutations['Variant_score'].to_list())
        self.df_single_mutations['Variant_class'] = np.select([mut_scores >= self.cluster_means[1], (mut_scores < self.cluster_means[1]) & (mut_scores >= self.thresh)], [1, 2], default=3)
        self.df_single_mutations['Variant_class'] = self.df_single_mutations['Variant_class'].replace({1: 'neutral', 2: 'mild', 3: 'impactful'})
        residue_counts = self.df_single_mutations.groupby('Residue')['Variant_class'].value_counts(normalize=True).unstack(fill_value=0)
        residues_of_interest = residue_counts[(residue_counts['mild'] + residue_counts['impactful']) >= 0.5].index.tolist()
        self.df_single_mutations['Residue_class'] = 'tolerant'
        self.df_single_mutations.loc[self.df_single_mutations['Residue'].isin(residues_of_interest), 'Residue_class'] = 'sensitive'
        self.df_single_mutations.to_csv(self.GEMMEmat_location, index=False)



    #%%
    def plotGMM(self):
    
        '''plot GMM distribution for a given protein and its mutations'''
        #sns.set(style="ticks")
        #sns.set_context("notebook")
        
        bin_width = 0.25
        data_range = np.max(self.X_train) - np.min(self.X_train)
        num_bins = int(np.ceil(data_range / bin_width))
        bin_edges = np.linspace(np.min(self.X_train), np.max(self.X_train), num_bins + 1)
        nx_pmf, xbins = np.histogram(self.X_train, bins=bin_edges,density=False)
        nx_pmf = nx_pmf/len(self.df_single_mutationsConf)
        widthBar = xbins[1] - xbins[0]
        stdD = np.sqrt(np.array(self.model_GMM.covariances_).flatten()[self.sorted_indices])

        fig, ax = plt.subplots(figsize=(8, 5)) 
        #fig, ax=plt.figure(figsize=(8, 5))
        plt.bar(xbins[:-1], nx_pmf, width=widthBar, align='edge', color='xkcd:grey',edgecolor='white', alpha=0.5, label='Variant scores')
        xbins[-1] = xbins[-1] + 0.00001
        widthBar = xbins[1] - xbins[0]
        cluster_colors = ['xkcd:red', 'xkcd:purple','xkcd:blue']
        clusters_height = pd.DataFrame({'xbins':xbins[:-1]})
        self.component_shareOrdered = np.hstack((self.X_train, self.component_shareOrdered))
        # Plot histograms for each cluster
        for cluster, color in zip(range(1, self.component_shareOrdered.shape[1]),cluster_colors):  # Iterate over clusters, not probabilities
            l_height =[np.mean(self.component_shareOrdered[(self.component_shareOrdered[:, 0] >= xbins[i-1]) &
                                                        (self.component_shareOrdered[:, 0] < xbins[i]), cluster])
                        for i in range(1, len(xbins))] * nx_pmf
            clusters_height[cluster-1]=l_height
            label_custom=r"$\mu={:.2f} \ ; \ \sigma={:.2f}$".format(self.cluster_means[cluster-1],stdD[cluster-1])
            plt.bar(xbins[:-1], l_height, color=color, alpha=0.27, label=label_custom,edgecolor='white',  width=widthBar, align='edge')   

        plt.axvline(x=self.thresh, color='#6f0043', linestyle='-', linewidth=3.5) 
        plt.axvline(x=self.cluster_means[1], color='#ffcfc4', linestyle='-', linewidth=3.5)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        #plt.xticks(fontsize=27, weight='bold')
        #plt.yticks(fontsize=27)
        plt.ylabel('Density', fontsize=13, labelpad=10)
        plt.xlabel('Variant score', fontsize=13, labelpad=10)
        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=10)
        #plt.legend(loc = 'upper left', fontsize=23)
        plt.tight_layout()
        plt.savefig(f'{self.outputDirPickle}/{self.prot_name}_GMM_model.jpg', dpi=500, format='jpg')
        plt.close()

if __name__ == '__main__':
    proteocast = sys.argv[1]
    protein = '_'.join(proteocast.split('/')[-1].split('_')[:-1])
    plotGMM = sys.argv[2].lower() == 'true'

    
    GMM(GEMMEmat_location=proteocast, prot_name=protein, nb_clusters=3, outPickle='/'.join(proteocast.split('/')[:-1])+'/GMM/', plot=plotGMM)
