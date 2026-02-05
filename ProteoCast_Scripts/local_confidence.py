### calculate the local confidence for the whole proteome using the GEMME predictions
### weighted smoothing with a gaussian kernel (convolution) for intervals

import numpy as np 
from utils import read_GEMME_mat,removeWTvalue, alph
import pandas as pd
from Bio import SeqIO
import warnings
from pathlib import Path
import argparse
import os
import sys

warnings.filterwarnings("ignore", category=FutureWarning)

WORK_DIR = Path.cwd()

def convert_to_intervals(lst):
    if lst==[]:
        return []
    intervals = []
    start = lst[0]
    for i in range(1, len(lst)):
        if lst[i] != lst[i - 1] + 1:
            end = lst[i - 1]
            if start != end:  # Exclude intervals where start == end
                intervals.append((start, end))
            else:
                intervals.append(lst[i-1])
            start = lst[i]
    end = lst[-1]
    if start != end:  # Exclude intervals where start == end
        intervals.append((start, end))
    else:
        intervals.append(end)
    return intervals



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="GEMME results confidence analysis.")
    parser.add_argument("--gemme_dir", required=True, help="Path to GEMME predictions.")
    parser.add_argument("--msa_name", required=True, help="MSA name.")
    args = parser.parse_args()

    if not args.gemme_dir:
        print("Error: The '--gemme_dir' argument is missing.")
        parser.print_help()
    elif not args.msa_name:
        print("Error: The '--msa_name' argument is missing")
        parser.print_help()
    else:
        protein_id = args.gemme_dir.split('/')[-1]
        msa_name = args.msa_name

        pad_value = 1  # Value to pad the sequence with
        weights = [2., 3., 4., 3., 2.]
        std_thresh = 0.15; cons_thresh = 0.2; score_thresh = 0.1

        try:
            mat = read_GEMME_mat(args.gemme_dir, protein_id, subdir=False)
            mat_conservation = pd.read_csv(os.path.join(args.gemme_dir, f'{protein_id}_conservation.txt'),sep=' ')        ##pos 250 -> 249th element 
            mat_Epi = read_GEMME_mat(args.gemme_dir, protein_id, mat_name='_pred_evolEpi.txt',subdir=False)      
        except Exception as e:
            print(f"Error: {e}")
            sys.exit(1)

        if mat.min() == 0.0:
            print("Matrix contains only zero values.")
            sys.exit(1) 

        seq = list(SeqIO.parse(args.gemme_dir+'/'+protein_id+'.fasta', 'fasta'))[0].seq
        MSA = list(SeqIO.parse('/'.join(args.gemme_dir.split('/')[:-1])+f'/ali{msa_name}.fasta', 'fasta'))

        df_single_mutations = pd.DataFrame(columns=['Mutation', 'Variant_score'])
        df_single_mutations['Mutation'] = [WT+str(i+1)+MUT for i, WT in enumerate(seq)for MUT in alph[:20]]
        df_single_mutations['Variant_score'] = mat.flatten(order='F')
        df_single_mutations['Residue'] = np.repeat(np.arange(1, len(seq)+1), 20)
       
        if len(MSA)<200:
            df_single_mutations['LocalConfidence'] = False
        else:
            df_single_mutations['LocalConfidence'] = True
            

        ## std in a column
        new_mat = removeWTvalue(mat);std_values = np.std(new_mat, axis=0)
        f_obs = 1-(np.sum(mat_Epi == 1, axis=0)/19); coverage = 1-mat_conservation.loc['gap']
        conservation= mat_conservation.loc['trace']

        ConfidanceLevel = (f_obs)*(coverage)
        
        residues_LowConf = ConfidanceLevel[ConfidanceLevel < score_thresh].index.astype(int).tolist()
        residues_LowConf = [x for x in residues_LowConf if (std_values[x-1]< std_thresh and conservation.loc[str(x)] < cons_thresh)]
        if residues_LowConf!=[]:
            mat_confidence = [1 if i in residues_LowConf else 0 for i in range(1, mat.shape[1] + 1)]
            ## a paddeed version for the extremeties
            mat_confidence = np.pad(mat_confidence, (2, 2), mode='constant', constant_values=pad_value)

            ## a paddeed version for the extremities
            seq_convolve = np.convolve(mat_confidence, weights, 'same')/sum(weights)
            ConfidanceLevel_Convolution = [int(x+1) for x in np.where(np.array(seq_convolve[2:-2])>=0.5)[0]]
            for pos in ConfidanceLevel_Convolution:
                df_single_mutations.loc[df_single_mutations['Residue']==pos, 'LocalConfidence'] = False


        df_single_mutations.to_csv('/'.join(args.gemme_dir.split('/')[:-1]) + '/' + protein_id + '_ProteoCast.csv', index=False)


