import pandas as pd 
import numpy as np
import gemmi as gm
from Bio import SeqIO
from prody import *
import matplotlib.pyplot as plt
import seaborn as sns
import pickle

#################    Data  ####################

alph = ["a","c","d","e","f","g","h","i","k","l","m","n","p","q","r","s","t","v","w","y","-",'x']
alph = [i.upper() for i in alph]
amino_acids = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
df_Dmel_recap = pd.read_csv('../csv/Dmel6.44PredictionsRecap.csv', index_col=0)
AF_dir = '../AF_structures/'



#################   GEMME data manipulation ####################
#%%
def read_GEMME_mat(dir_name, FBpp_ID, mat_name='_normPred_evolCombi.txt', subdir=True):
    """
    Read GEMME matrix; Combi by default or Ind and Epi if given.
    
    Parameters:
    dir_name (str): Directory name where the GEMME matrix is located.
    FBpp_ID (str): Protein name.
    mat_name (str): Matrix file name. Default is '_normPred_evolCombi.txt'.
    subdir (bool): If True, looks for the matrix in a subdirectory named after FBpp_ID. Default is True.
    
    Returns:
    numpy.ndarray: The GEMME matrix.
    """
    if subdir:
        f_path = f"{dir_name}/{FBpp_ID}/{FBpp_ID}{mat_name}"
    else:
        f_path = f"{dir_name}/{FBpp_ID}{mat_name}"
    
    mat = np.loadtxt(f_path, dtype='str', skiprows=1)
    mat = np.delete(mat, 0, axis=1)  # Delete the letters at the beginning of each matrix's row
    mat[mat == 'NA'] = 0.0
    mat = mat.astype(float)

    return mat
#%%
def removeWTvalue(data):
    """
    Removes the first occurrence of 0.0 from each column in the input data matrix.

    Parameters:
    data (numpy.ndarray): A 2D numpy array from which the first occurrence of 0.0 
    data (numpy.ndarray): A 2D numpy array from which the first occurrence of 0.0 
                          in each column will be removed.

    Returns:
    numpy.ndarray: A new 2D numpy array with the same number of columns as the input 
                   data, but with one less row, where the first occurrence of 0.0 
                   in each column has been removed.
    """
    new_mat = np.zeros((19, data.shape[1]))
    # Remove only one occurrence of 0.0 from each column
    for col_index in range(data.shape[1]):
        col = list(data[:, col_index])
        col.remove(0.0)
        new_mat[:, col_index] = col
    return new_mat

#%%
def ResidueSensitivityGEMME(FBpp_ID, l_residue=[]):
    """
    Analyze the sensitivity of residues in a given protein based on GEMME predictions.
    This function classifies residues of a protein into different sensitivity categories
    using thresholds for impactful and uncertain classifications. It then identifies 
    residues that are highly sensitive.
    Parameters:
    FBpp_ID (str): The FlyBase protein identifier for the protein of interest.
    l_residue (list, optional): A list of residue positions to check for high sensitivity. 
                                If provided, the function returns a list indicating whether 
                                each residue in l_residue is highly sensitive. Defaults to an empty list.
    Returns:
        If l_residue is provided, returns a list of integers where each 
        integer indicates whether the corresponding residue in l_residue 
        is highly sensitive (1) or not (0). If l_residue is not provided, 
        returns a numpy array of residue positions that are highly sensitive.
    """

    # Get the thresholds for impactful and uncertain classifications
    impactful_threshA = df_Dmel_recap.loc[FBpp_ID, 'GMM3_impactful']
    uncertain_threshA = df_Dmel_recap.loc[FBpp_ID, 'GMM3_uncertain']
    
    # Read GEMME predictions for the given protein
    matCombiA = read_GEMME_mat('../CombiPred/', FBpp_ID, subdir=False)
    mu_mat = -np.mean(matCombiA, axis=0)
    mu_mat = (mu_mat - np.min(mu_mat))/(np.max(mu_mat) - np.min(mu_mat))
    # Classify mutation scores based on thresholds
    classificationsA = np.select(
        [matCombiA >= uncertain_threshA, (matCombiA < uncertain_threshA) & (matCombiA >= impactful_threshA)], [1, 2], default=3)

    # Count the number of residues classified as uncertain or impactful in each column
    count_2_or_3 = np.sum((classificationsA == 2) | (classificationsA == 3), axis=0)

    # Identify residues that are highly sensitive (classified as 2 or 3 in at least 10 columns)
    highly_sensitive = np.where(count_2_or_3 >= 10)[0] + 1
    
    # If a list of residues is provided, return a list indicating whether each residue is highly sensitive
    if l_residue:
        return np.array([int(residue in highly_sensitive) for residue in l_residue])
    # Otherwise, return the list of highly sensitive residues
    else:
        return highly_sensitive, mu_mat
#%%
############################ 3D structure manipulation ############################

#%%
def bfactorsPDB(modelPDB):
    df_bfactors = pd.DataFrame(columns=['bfactor', 'aa'])
    for resi in modelPDB:
        for atom in resi:
            df_bfactors.loc[resi.seqid.num, 'bfactor'] = round(atom.b_iso,3)
            df_bfactors.loc[resi.seqid.num, 'aa'] = amino_acids[resi.name]
            break
    return df_bfactors

#%%
def PDBtoSeq(uniprot):
    model = gm.read_pdb(AF_dir+f'AF-{uniprot}-F1-model_v4.pdb')[0][0]
    sequence=''.join(amino_acids[residue.name] for residue in model)

    return sequence

# %%
def Change_Bfactors(pdb_path, gemme_path, new_pdb_path, column):
    '''
    column (str): 'max' or 'mean' 
    change the structure's bfactors to mean or max score per residue
    '''

    # Read a PDB file
    doc = gm.read_pdb(pdb_path)
    gemme_mat = read_GEMME_mat(gemme_path,'rank')

    if column=='max':
        bfactors = gemme_mat.max(axis=0)
    elif column=='mean':
        bfactors = gemme_mat.mean(axis=0)
    
    for chain in doc[0]:
        cpt=0
        for res in chain:
            for at in res:
                at.b_iso = bfactors[cpt]
            cpt+=1
    # # Write a PDB file
    doc.write_pdb(new_pdb_path)
#%%
############################    Plots   ############################
#%%
def plotGMM(FBpp, l_mut=[], SNPscores=[], SNPset = [], path=None):
    
    '''plot GMM distribution for a given protein and its mutations'''
    sns.set(style="ticks")
    sns.set_context("notebook")

    GMM_model=pickle.load(open('../GMM3_Confidence'+f'/{FBpp}_GMM_model.pickle', 'rb'))
    Single_mut=pd.read_csv('../GMM3_Confidence_Labels'f'/{FBpp}_labels.csv', sep=',')
    # Extract GEMME raw values from Single_mut DataFrame
    values = np.array(Single_mut['gemme_score'])

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
    nx_pmf = nx_pmf/len(Single_mut)
    widthBar = xbins[1] - xbins[0]

    fig, ax = plt.subplots(figsize=(16, 10)) ##18, 14
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
            dc_colors = {1:'r', 0.5:'m', 0:'b'}
            dc_colors_set = {'Lethal':'r', 'DGRP':'b'}
            #SNPscores = SNPscores[sorted_mut]; l_mut = l_mut[sorted_mut];l_GMMclass = l_GMMclass[sorted_mut]
            for i in range(len(l_mut)):
                #color = dc_colors[SNPclasses[sorted_mut[i]]]
                color = dc_colors_set[SNPset[sorted_mut[i]]]
 
                if SNPscores[sorted_mut[i]]<df_Dmel_recap.loc[FBpp, 'GMM3_impactful']:
                    y = max(nx_pmf)/2- 0.003*i
                else:
                    y = max(nx_pmf)- 0.003*i 
                #if l_mut[sorted_mut[i]]=='P184S':
                #    y = max(nx_pmf)- 0.01
                y_line = np.arange(0, y, 0.001)
                # Plot dots and lines for mutations
                plt.plot(SNPscores[sorted_mut[i]], y, f'o{color}--', alpha=0.7, markersize=11)
                plt.plot([SNPscores[sorted_mut[i]]] * len(y_line), y_line, 'k:', alpha=0.3, linewidth=3)
                # Annotate mutations
                plt.annotate(f'{l_mut[sorted_mut[i]]}', (SNPscores[sorted_mut[i]] + 0.045, y), fontsize=27, color='black',alpha=0.55, weight='normal',
                                rotation=30)#dimgrey






    plt.axvline(x=df_Dmel_recap.loc[FBpp, 'GMM3_impactful'], color='#6f0043', linestyle='-', linewidth=4.5) #label = ' = - 2.25',
    plt.axvline(x=df_Dmel_recap.loc[FBpp, 'GMM3_uncertain'], color='#ffcfc4', linestyle='-', linewidth=4.5) #label = ' threshold = - 3.48',

    # Hide all spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.xticks(fontsize=27, weight='bold')
    plt.yticks(fontsize=27)
    plt.ylabel('Density', fontsize=37, labelpad=25)
    plt.xlabel('GEMME score', fontsize=37, labelpad=29)
    ## alternative
    ax.tick_params(axis='x',pad=10)

    plt.legend(loc = 'upper left', fontsize=23)
    plt.tight_layout()
    #plt.title(r'GMM distribution for $\mathit{ yorkie}$ prEVEotein ', fontsize=20 )
    if path:
        plt.savefig(path, dpi=500, format='jpg')
    #plt.show()
    plt.close()


#%%
## PDB structure retreival (3D not included in the batch )
# import requests

#l_noPDBstruct = []
#for uniprot in df_Uniprot_FBpp.loc[~df_Uniprot_FBpp['UniProt'].isin(l_3Dstructure), 'UniProt'].unique():
#
#    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v4.pdb"
#
#    # Fetch and save the PDB file
#    response = requests.get(url)
#
#    if response.status_code == 200:
#        with open('../AF_structures_outOfBulk/'+f"AF-{uniprot}-F1-model_v4.pdb", "wb") as file:
#            file.write(response.content)
#    else:
#        l_noPDBstruct.append(uniprot)
