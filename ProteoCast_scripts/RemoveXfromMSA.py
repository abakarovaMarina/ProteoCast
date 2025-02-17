import os
from Bio import SeqIO
from pathlib import Path
import argparse
import shutil

WORK_DIR = Path.cwd()

def most_present(MSA, index):
    l_aa = [str(occ.seq)[index] for occ in MSA[1:]]
    most_common = max(set(l_aa), key=l_aa.count)
    return most_common

def rewriteFasta(pp_MSA, l_ind, f_out):
    seq_ppQ = pp_MSA[0].seq

    for ind in l_ind:
        aa_mostPresent = most_present(pp_MSA, ind)
        if aa_mostPresent == '-':
            for record in pp_MSA:
                record.seq = record.seq[:ind] + record.seq[ind+1:]
        else:
            seq_ppQ = seq_ppQ[:ind]+aa_mostPresent+seq_ppQ[ind+1:]
            pp_MSA[0].seq = seq_ppQ
    
    with open(f_out, "w") as f:
        SeqIO.write(pp_MSA, f, "fasta")
    return 1

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Remove undefined amino acids from the MSA.")
    parser.add_argument("--msa", required=True, help="Path to the input msa fasta file.")
    args = parser.parse_args()

    if not args.msa:
        print("Error: The '--msa' argument is missing.")
        parser.print_help()
    else:
        protein_id = args.msa.split('/')[-1][3:-6]
        msa_path = os.path.join(WORK_DIR, args.msa)
        
        if os.path.exists(msa_path):
            pp_MSA = list(SeqIO.parse(args.msa, "fasta"))
            seq_ppQ = pp_MSA[0].seq
            l_ind = [i for i in range(len(seq_ppQ)) if seq_ppQ[i] == 'X']
            
            if l_ind:
                print('there is X')
                f_out = '/'.join(args.msa.split('/')[:-1]) + '/ali' + protein_id + 'X.fasta'
                shutil.move(args.msa, f_out)
                rewriteFasta(pp_MSA, l_ind, args.msa)
        else:
            print('The provided fasta file does not exist')