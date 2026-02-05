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

#df_Dmel_recap = pd.read_csv('../csv/Dmel6.44PredictionsRecap.csv', index_col=0)
AF_dir = 'AF_structures/'



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
        return highly_sensitive
#%%
############################ 3D structure manipulation ############################

#%%
def bfactorsPDB(modelPDB):
    df_bfactors = pd.DataFrame(columns=['pLDDT', 'AA'])
    for resi in modelPDB:
        for atom in resi:
            df_bfactors.loc[resi.seqid.num, 'pLDDT'] = round(atom.b_iso,3)
            df_bfactors.loc[resi.seqid.num, 'AA'] = amino_acids[resi.name]
            break
    return df_bfactors

#%%
def PDBtoSeq(uniprot):
    model = gm.read_pdb(AF_dir+f'AF-{uniprot}-F1-model_v4.pdb.gz')[0][0]
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
    gemme_mat = read_GEMME_mat(gemme_path)

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

