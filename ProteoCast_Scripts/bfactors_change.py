# ### JUPYTER NOTEBOOK TO MAP GEMME SCORES ON THE PROTEIN 3D STRUCTURE 

# %%
import gemmi as gm
import numpy as np
import argparse
import pandas as pd


def Change_Bfactors_MappedChain(pdb_path, chain_id, bfactor_map, new_pdb_path, default_bfactor=np.nan):
    '''
    Updates B-factors on a specific chain based on a residue → B-factor mapping.
    Unmapped residues are left as default_bfactor (NaN or 0.0).

    Parameters:
    - pdb_path: Path to input PDB
    - chain_id: e.g., 'A'
    - bfactor_map: dict or Series with {residue_number: bfactor}
    - new_pdb_path: Path to output modified PDB
    - default_bfactor: Value to use for unmapped residues (NaN will not write b_iso)
    '''
    print("output pdb:", new_pdb_path)
    doc = gm.read_pdb(pdb_path)
    model = doc[0]
    for chain in model:
        for residue in chain:
            resid = residue.seqid.num

            if chain.name == chain_id and resid in bfactor_map:
                bfactor = bfactor_map[resid]
            else:
                bfactor = 0.0  # default for all other residues

            for atom in residue:
                atom.b_iso = float(bfactor)

    doc.write_pdb(new_pdb_path)



if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Load a PDB structure and add custom B-factors.")
    parser.add_argument("--pdb", required=True, help="Path to the pdb.")
    parser.add_argument("--mapping", required=True, help="Path to the mapping file.")
    parser.add_argument("--chain", required=True, help="Chain identifier.")
    parser.add_argument("--rsa", required=True, help="Path to the rsa file.")
    parser.add_argument("--proteocast", required=True, help="Path to the proteocast file.")

    args = parser.parse_args()

    if not args.pdb:
        print("Error: The '--pdb' argument is missing.")
        parser.print_help()
    elif not args.mapping:
        print("Error: The '--mapping' argument is missing.")
        parser.print_help()
    elif not args.chain:
        print("Error: The '--chain' argument is missing.")
        parser.print_help()
    elif not args.rsa:
        print("Error: The '--rsa' argument is missing.")
        parser.print_help()
    elif not args.proteocast:
        print("Error: The '--proteocast' argument is missing.")
        parser.print_help()
    else:
        df_proteocast = pd.read_csv(args.proteocast, sep=',')
        df_mapping = pd.read_csv(args.mapping, sep=',')
        df_rsa = pd.read_csv(args.rsa, sep=',')
        print(args.pdb, args.chain)

        ## add 'aligned' column to rsa_values.csv dataframe for visualization purposes 
        # initialize all residues as unaligned 
        df_rsa['aligned'] = False
        #set aligned = True for residues in the given chain and in the mapping
        mask = (df_rsa['chain'] == args.chain) & (df_rsa['residue'].isin(df_mapping['pdb_pos']))
        df_rsa.loc[mask, 'aligned'] = True
        df_rsa.to_csv(args.rsa, index=False)

        ## RSA*Variant_score mapping
        # 1. Filter RSA for the relevant chain
        df_rsa_chain = df_rsa[df_rsa['chain'] == args.chain]
        # 2. Merge ProteoCast with mapping on FASTA position
        df_merged = pd.merge(df_proteocast, df_mapping, left_on="Residue", right_on="fasta_pos", how="inner")
        # 3. Merge with RSA using PDB residue index
        df_merged = pd.merge(df_merged, df_rsa_chain, left_on="pdb_pos", right_on="residue", how="inner")
        # 4. Compute new column
        df_merged["RSA*Variant_score"] = ((1.0 - df_merged["rsa"].clip(upper=100) / 100.0)* df_merged["Variant_score"])

        # 5. (optional) push back into df_proteocast if you want to align on Residue
        df_proteocast = pd.merge(df_proteocast, df_merged[["Mutation", "pdb_pos", "RSA*Variant_score"]], on="Mutation", how="left")
        df_proteocast['aligned'] = df_proteocast["RSA*Variant_score"].notna()
        #df_proteocast["RSA*Variant_score"] = df_proteocast["RSA*Variant_score"].fillna(df_proteocast["Variant_score"])

        # 6. Save the merged DataFrame
        df_proteocast.to_csv(args.proteocast, index=False)

        # Group and average over all variants at the same PDB residue
        bfactor_grouped = df_proteocast.groupby("pdb_pos").agg({"RSA*Variant_score": "mean","Variant_score": "mean","Residue_class": "first"})

        # Normalize *after* averaging
        bfactor_grouped["RSA*Variant_score_norm"] = -bfactor_grouped["RSA*Variant_score"]
        min_rsa = bfactor_grouped["RSA*Variant_score_norm"].min()
        max_rsa = bfactor_grouped["RSA*Variant_score_norm"].max()
        bfactor_grouped["RSA*Variant_score_norm"] = ((bfactor_grouped["RSA*Variant_score_norm"] - min_rsa) / (max_rsa - min_rsa) * 100)

        bfactor_grouped["Variant_score_norm"] = -bfactor_grouped["Variant_score"]
        min_var = bfactor_grouped["Variant_score_norm"].min()
        max_var = bfactor_grouped["Variant_score_norm"].max()
        bfactor_grouped["Variant_score_norm"] = (
            (bfactor_grouped["Variant_score_norm"] - min_var) / (max_var - min_var) * 100)

        # Map to dictionaries
        bfactor_ResClass = bfactor_grouped["Residue_class"].map({"sensitive": 100, "tolerant": 0}).to_dict()
        bfactor_rsa = bfactor_grouped["RSA*Variant_score_norm"].to_dict()
        bfactor_sensitivity = bfactor_grouped["Variant_score_norm"].to_dict()
      
        out_pdbResClass = 'job/11.'+args.pdb[:-4].split('/')[1] + "_ResClass.pdb"
        out_pdbSensitivity = 'job/12.'+args.pdb[:-4].split('/')[1] + "_Sensitivity.pdb"
        out_pdbrsa = 'job/13.'+args.pdb[:-4].split('/')[1] + "_ResRSA.pdb"
        
        # Apply it
        Change_Bfactors_MappedChain(args.pdb, args.chain, bfactor_rsa, out_pdbrsa)
        Change_Bfactors_MappedChain(args.pdb, args.chain, bfactor_ResClass, out_pdbResClass)
        Change_Bfactors_MappedChain(args.pdb, args.chain, bfactor_sensitivity, out_pdbSensitivity)

    

        
