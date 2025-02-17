# %% [markdown]
# ### JUPYTER NOTEBOOK TO MAP GEMME SCORES ON THE PROTEIN 3D STRUCTURE 

# %%
import gemmi as gm
import numpy as np
import argparse
import pandas as pd


# %%
def Change_Bfactors(pdb_path, bfactors, new_pdb_path):
    '''
    column (str): 'max' or 'mean' 
    change the structure's bfactors to mean or max score per residue
    '''
    print(max(bfactors), min(bfactors))
    # Read a PDB file
    doc = gm.read_pdb(pdb_path)
    for chain in doc[0]:
        cpt=0
        for res in chain:
            for at in res:
                at.b_iso = bfactors[cpt]
            cpt+=1
    # # Write a PDB file
    doc.write_pdb(new_pdb_path)

#%%
def Change_Bfactors_Sensitivity(pdb_path, l_sensitiveRes, new_pdb_path):
    '''
    column (str): 'max' or 'mean' 
    change the structure's bfactors to mean or max score per residue
    '''

    # Read a PDB file
    doc = gm.read_pdb(pdb_path)

    for chain in doc[0]:
        cpt=0
        for res in chain:
            res_num = res.seqid.num
            if res_num in l_sensitiveRes:
                bfactor=100
            else:
                bfactor=0
            for at in res:
                at.b_iso = bfactor
            cpt+=1
    # # Write a PDB file
    doc.write_pdb(new_pdb_path)
#%% 

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Load a PDB structure and add custom B-factors.")
    parser.add_argument("--pdb", required=True, help="Path to the pdb.")
    parser.add_argument("--proteocast", required=True, help="Path to the proteocast file.") 
    args = parser.parse_args()

    if not args.pdb:
        print("Error: The '--pdb' argument is missing.")
        parser.print_help()
    elif not args.proteocast:
        print("Error: The '--proteocast' argument is missing.")
        parser.print_help()
    else:
        df_proteocast = pd.read_csv(args.proteocast, sep=',')
        print(args.pdb)
        l_sensitive = df_proteocast.loc[df_proteocast['Residue_class']=='sensitive', 'Residue'].unique().tolist()
        structure = args.pdb
        out_pdbResClass = args.pdb[:-4] + "_GEMMEResClass.pdb"
        Change_Bfactors_Sensitivity(args.pdb,l_sensitive, out_pdbResClass)

        mu_mat = - df_proteocast.groupby('Residue')['GEMME_score'].mean()
        bfactors = list((mu_mat - np.min(mu_mat))/(np.max(mu_mat) - np.min(mu_mat))*100)
        out_pdbSensitivity = args.pdb[:-4] + "_GEMMESensitivity.pdb"
        Change_Bfactors(args.pdb, bfactors, out_pdbSensitivity)
