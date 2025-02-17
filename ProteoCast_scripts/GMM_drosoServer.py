import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import random
import matplotlib.pyplot as plt
import sys
import warnings
from sklearn.exceptions import InconsistentVersionWarning


warnings.filterwarnings("ignore", category=InconsistentVersionWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

os.chdir('/data/')

def plotGMM(GMM_model, df_thresh, df_proteocast, l_mut=[], SNPscores=[], SNPset = [], path=None):
   
    '''plot GMM distribution for a given protein and its mutations'''
    #sns.set(style="ticks")
    #sns.set_context("notebook")
   
    # Extract GEMME raw values from Single_mut DataFrame
    values = np.array(df_proteocast['GEMME_score'])

    # Calculate component shares (probabilities) using GMM for the given data values
    component_share = GMM_model.predict_proba(values.reshape(-1, 1))  # Shape: (#values, #clusters)
    # Get the cluster means and sort the indices
    cluster_means = GMM_model.means_.flatten(); sorted_indices = np.argsort(cluster_means)
    cluster_means = cluster_means[sorted_indices]; stdD = np.sqrt(np.array(GMM_model.covariances_).flatten()[sorted_indices])
    # Reorder component_share based on sorted indices and add the values column
    component_shareOrdered = component_share[:, sorted_indices]
    component_shareOrdered = np.hstack((values.reshape(-1, 1), component_shareOrdered))
    bin_width = 0.25
    # Determine the desired bin width
    # Calculate the number of bins dynamically based on the data range and desired bin width
    data_range = np.max(values) - np.min(values)
    num_bins = int(np.ceil(data_range / bin_width))
    bin_edges = np.linspace(np.min(values), np.max(values), num_bins + 1)
    nx_pmf, xbins = np.histogram(values, bins=bin_edges,density=False)
    nx_pmf = nx_pmf/len(df_proteocast)
    widthBar = xbins[1] - xbins[0]

    fig, ax = plt.subplots(figsize=(8, 5)) ##18, 14
    plt.bar(xbins[:-1], nx_pmf, width=widthBar, align='edge', color='xkcd:grey',edgecolor='white', alpha=0.5, label='GEMME scores')
    #plt.hist(Single_mut, color='xkcd:grey', bins=bin_edges, histtype='stepfilled', alpha=0.5, density=True)#density= True
    # Adjust the last bin edge
    xbins[-1] = xbins[-1] + 0.00001
    # Define colors and width for bars
    widthBar = xbins[1] - xbins[0]
    # Define colors and width for bars
    set_colors = ['xkcd:purple', 'xkcd:green', 'xkcd:grey']                                             ## colors for 3,4,5 GMM clusters
    # Get colors for the clusters (skip the first and last colors)
    cluster_colors = ['xkcd:red']+[set_colors[i] for i in range(len(cluster_means- 1)-2)]+['xkcd:blue']
    clusters_height = pd.DataFrame({'xbins':xbins[:-1]})

    dc_colors_set = {'Lethal':'r', 'DGRP':'b', 'DEST2':'b', 'DGRP+DEST2':'b', 'DEST2+DGRP':'b'}
    marker_shapes = {'Lethal': 'o', 'DGRP': 's', 'DEST2': '^', 'DEST2+DGRP': 'X', 'DGRP+DEST2': 'X'}
    # Plot histograms for each cluster
    for cluster, color in zip(range(1, component_shareOrdered.shape[1]),cluster_colors):  # Iterate over clusters, not probabilities
        # Calculate heights for the bars
        l_height =[np.mean(component_shareOrdered[(component_shareOrdered[:, 0] >= xbins[i-1]) &
                                                    (component_shareOrdered[:, 0] < xbins[i]), cluster])
                    for i in range(1, len(xbins))] * nx_pmf
        clusters_height[cluster-1]=l_height
        # Plot the bars uncomment if want to plot
        label_custom=r"$\mu={:.2f} \ ; \ \sigma={:.2f}$".format(cluster_means[cluster-1],stdD[cluster-1])
        plt.bar(xbins[:-1], l_height, color=color, alpha=0.27, label=label_custom,edgecolor='white',  width=widthBar, align='edge')
        
        # Iterate over mutations and add dots and annotations
        if SNPscores:
            sorted_mut = np.argsort(SNPscores)
            
            cpt=0; flag=True
            for i in range(len(l_mut)):
                #color = dc_colors[SNPclasses[sorted_mut[i]]]
                color = dc_colors_set[SNPset[sorted_mut[i]]]
                x_pos = SNPscores[sorted_mut[i]]
        
                if  x_pos < np.min(values) / 2:
                    y = max(nx_pmf) / 2 - 0.003 * i
                elif x_pos > np.min(values) / 4:
                    if (x_pos>-1) and flag:
                        cpt=0;flag=False
                    y = max(nx_pmf) - 0.0045 * cpt
                    
                    cpt += 1
                else:
                    y = max(nx_pmf) - 0.003 * i

                y_line = np.arange(0, y, 0.001)
                # Plot dots and lines for mutations
           
                plt.plot(x_pos, y, marker=marker_shapes[SNPset[sorted_mut[i]]], color=color, zorder=3, alpha=0.5, markersize=4, markeredgewidth=1)
                plt.plot([x_pos] * len(y_line), y_line, 'k:', alpha=0.2, linewidth=1, zorder=2)
                # Annotate mutations
                plt.annotate(f'{l_mut[sorted_mut[i]]}', (x_pos + 0.045, y), fontsize=9, color='black', alpha=0.55, weight='normal',
                             rotation=45, zorder=4, bbox=dict(boxstyle="round,pad=0.3", edgecolor="none", facecolor="white", alpha=0.4))


    plt.axvline(x=df_thresh.loc['GMM_impactful'].item(), color='#6f0043', linestyle='-', linewidth=3, zorder=1) #label = ' = - 2.25',
    plt.axvline(x=df_thresh.loc['GMM_uncertain'].item(), color='#ffcfc4', linestyle='-', linewidth=3, zorder=1) #label = ' threshold = - 3.48',

    # Hide all spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.ylabel('Density', fontsize=13, labelpad=10)
    plt.xlabel('GEMME score', fontsize=13, labelpad=10)
    ax.tick_params(axis='x', labelsize=10)
    ax.tick_params(axis='y', labelsize=10)
    ## alternative
    #ax.tick_params(axis='x',pad=10)
    plt.tight_layout()
    if path:
        plt.savefig(path, dpi=500, format='jpg')
    plt.clf()
    plt.close()

if __name__ == '__main__':
    proteocast = sys.argv[1]
    df_proteocast = pd.read_csv(proteocast, index_col=0)

    protein = proteocast.split('/')[-1].split('_')[0].split('.')[1]
    id = int(proteocast.split('/')[-2])

    df_thresh = pd.read_csv(f'Drosophila_ProteoCast/{id}/{protein}_GMMthresholds.csv', index_col=0)

    with open(f'Drosophila_ProteoCast/GMM3_Confidence/{protein}_GMM_model.pickle', 'rb') as f:
        GMM_model = pickle.load(f)

    if os.path.exists(f'Drosophila_ProteoCast/{id}/7.{protein}_SNPs.csv'):

        df_snps = pd.read_csv(f'Drosophila_ProteoCast/{id}/7.{protein}_SNPs.csv')
        df_snps = df_snps.loc[df_snps['Set_name']!='Hypomorphic'] 
        df_snps = df_snps.groupby(['Mutation', 'Residue']).agg({
                'GEMME_score': 'first',
                'Set_name': lambda x: '+'.join(sorted(set(x)))
            }).reset_index()
        df_snps = df_snps.loc[df_snps['Set_name'].isin(['Lethal', 'DGRP', 'DEST2','DGRP+DEST2', 'DEST2+DGRP'])]

        # Keep all Lethal mutations
        df_snps_lethal = df_snps[df_snps['Set_name'] == 'Lethal']
        # Take up to 15 DEST2+DGRP mutations
        df_snps_dest2_dgrp = df_snps[df_snps['Set_name'] == 'DEST2+DGRP'].head(15)

        # Calculate remaining space
        remaining_space = 20 - len(df_snps_dest2_dgrp)
        if remaining_space < 0:
            remaining_space = 0
        # Take remaining DEST2 mutations if space allows
        df_snps_dest2 = df_snps[df_snps['Set_name'] == 'DEST2'].head(remaining_space)

        # Update remaining space
        remaining_space -= len(df_snps_dest2)
        if remaining_space < 0:
            remaining_space = 0
        # Take remaining DGRP mutations if space allows
        df_snps_dgrp = df_snps[df_snps['Set_name'] == 'DGRP'].head(remaining_space)

        # Concatenate all selected mutations
        df_snps = pd.concat([df_snps_lethal, df_snps_dest2_dgrp, df_snps_dest2, df_snps_dgrp])


        #df_snps.drop_duplicates(subset=['Mutation', 'Residue'], keep='first', inplace=True)
        plotGMM(GMM_model,df_thresh, df_proteocast, protein, df_snps['Mutation'].tolist(), df_snps['GEMME_score'].tolist(),df_snps['Set_name'].tolist(), f'Drosophila_ProteoCast/{id}/6.{protein}_GMM.jpg')
    else: 
        plotGMM(GMM_model,df_thresh, df_proteocast, protein, path=f'Drosophila_ProteoCast/{id}/6.{protein}_GMM.jpg')

