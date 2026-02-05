import argparse
import pandas as pd
import numpy as np
import sys
import gemmi as gm


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="GEMME results confidence analysis.")
    parser.add_argument("--proteocast", required=True, help="Path to proteocast predictions file.")
    parser.add_argument("--mapping", help="Path to the mapping file.")
    parser.add_argument("--chain", help="Chain identifier.")
    parser.add_argument("--pdb", help="Path to pdb.")

    args = parser.parse_args()

    if not args.proteocast:
        print("Error: The '--proteocast' argument is missing.")
        parser.print_help()
    else:
        prot = '_'.join(args.proteocast.split('/')[-1].split('_')[:-1])
        try:
            df_proteocast = pd.read_csv(args.proteocast)   
        except Exception as e:
            print(f"Proteocast file was not found: {e}")
            sys.exit(1)

        
        mu_gemme = -df_proteocast.groupby('Residue')['Variant_score'].mean()
        alphafold = False
        if args.pdb:
            if 'AF' in args.pdb or 'alphafold' in args.pdb:
                alphafold = True
            else:
                with open(args.pdb, "r") as f:
                    for line in f:
                        if line.startswith("TITLE"):
                            header = line.strip()
                            if "ALPHAFOLD" in header:
                                alphafold = True
                                break

        if alphafold:                   
            df_mapping = pd.read_csv(args.mapping, sep=',')
            #every seq residue need to have a pLDDT value
            #mu_gemme = mu_gemme.loc[mu_gemme.index.isin(df_mapping['fasta_pos'].values)]
            GEMME_mean = (mu_gemme - np.min(mu_gemme))/(np.max(mu_gemme) - np.min(mu_gemme))

            model = gm.read_pdb(args.pdb)[0]
            chain = model[args.chain]  # Only your chain of interest

            bfactors = []
            for resi in chain:
                resid = resi.seqid.num
                if resid not in df_mapping['pdb_pos'].values:
                    continue  # skip unmapped residues

                # take only the first atom (any one works, e.g., CA)
                for atom in resi:
                    bfactors.append(round(atom.b_iso / 100, 3))
                    break
            #print('GEMME mean:', len(GEMME_mean), 'pLDDT:', len(bfactors))
            try:
                df = pd.DataFrame({'GEMME_mean': GEMME_mean, 'pLDDT': bfactors})
                df.to_csv('/'.join(args.proteocast.split('/')[:-1]) + '/' + prot + '_GEMME_pLDDT.csv') 
            except:
                print('different length of GEMME and pLDDT; GEMME', len(GEMME_mean), 'pLDDT', len(bfactors))
                GEMME_mean = (mu_gemme - np.min(mu_gemme))/(np.max(mu_gemme) - np.min(mu_gemme))
                df = pd.DataFrame({'GEMME_mean': GEMME_mean})
                df.to_csv('/'.join(args.proteocast.split('/')[:-1]) + '/' + prot + '_GEMME_pLDDT.csv')
        else:
            GEMME_mean = (mu_gemme - np.min(mu_gemme))/(np.max(mu_gemme) - np.min(mu_gemme))
            df = pd.DataFrame({'GEMME_mean': GEMME_mean})
            df.to_csv('/'.join(args.proteocast.split('/')[:-1]) + '/' + prot + '_GEMME_pLDDT.csv')
