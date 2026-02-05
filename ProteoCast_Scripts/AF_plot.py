"""PYTHON SCRIPT TO PRODUCE AN ALIGNMENT PLOT """

import numpy as np
import pandas as pd
from Bio import SeqIO
import os
from pathlib import Path
import gc
import matplotlib.pyplot as plt
import argparse
import sys

## amino acids alphabet
alph = ["a","c","d","e","f","g","h","i","k","l","m","n","p","q","r","s","t","v","w","y","-",'x', 'b', 'z', 'u', 'j']
alph = [i.upper() for i in alph]
d_aa_id = dict(zip(alph, range(1,len(alph)+1)))

## convert seqMSA to num MSA
def toNum(msa, d_aa_id):
    Ls = len(msa[0])
    N = len(msa)
    msaNum = np.array([[d_aa_id[x]for x in seq if x in alph] for seq in msa])
    return msaNum

## PLOT FUNCTION
def plot_msa_v2(MSA_file_path, d_aa_id, path, plotTitle, sort_lines=True, dpi=100):
    ## charge the MSA
    pp_MSA = [str(x.seq) for x in list(SeqIO.parse(MSA_file_path, "fasta"))]
    seqQuery = pp_MSA[0]
    seqQueryNum = [d_aa_id[x] for x in seqQuery]
    Ls = [len(seqQuery)]
    Ln = np.cumsum([0] + Ls)  ## cumulative length of sequences
    msaNum = toNum(pp_MSA, d_aa_id)
    gap = msaNum != 21  ## boolean indicating the presence of a gap ("-") in the MSA.
    qid = msaNum == seqQueryNum  ## boolean array where True indicates that the corresponding element in msa is equal to the query sequence

    gapid = np.stack([gap[:, Ln[i]:Ln[i + 1]].max(-1) for i in range(len(Ls))], -1)
    lines = []
    Nn = []

    for g in np.unique(gapid, axis=0):
        i = np.where((gapid == g).all(axis=-1))
        qid_ = qid[i]
        gap_ = gap[i]
        seqid = np.stack([qid_[:, Ln[i]:Ln[i + 1]].mean(-1) for i in range(len(Ls))], -1).sum(-1) / (g.sum(-1) + 1e-8)
        non_gaps = gap_.astype(float)
        non_gaps[non_gaps == 0] = np.nan
        if sort_lines:
            lines_ = non_gaps[seqid.argsort()] * seqid[seqid.argsort(), None]
        else:
            lines_ = non_gaps[::-1] * seqid[::-1, None]
        Nn.append(len(lines_))
        lines.append(lines_)

    Nn = np.cumsum(np.append(0, Nn))
    lines = np.concatenate(lines, 0)

    plt.figure(figsize=(8, 5), dpi=dpi)
    
    #plt.title(f"Sequence coverage", fontsize=15)
    plt.imshow(lines, interpolation='nearest', aspect='auto',
               cmap="rainbow_r", vmin=0, vmax=1, origin='lower',
               extent=(0, lines.shape[1], 0, lines.shape[0]))
    for i in Ln[1:-1]:
        plt.plot([i, i], [0, lines.shape[0]], color="black")
    for j in Nn[1:-1]:
        plt.plot([0, lines.shape[1]], [j, j], color="black")

    plt.plot((np.isnan(lines) == False).sum(0), color='black')
    plt.rcParams['savefig.facecolor'] = 'white'
    plt.xlim(0, lines.shape[1])
    plt.ylim(0, lines.shape[0])
    cbar = plt.colorbar(label="Sequence identity to query")
    cbar.set_label("Sequence identity to query", size=13,labelpad=7)
    plt.xlabel("Residue", fontsize=13, labelpad=10)
    plt.ylabel("Sequences", fontsize=13, labelpad=10)
    plt.tick_params(axis='y', labelsize=10)
    plt.tick_params(axis='x', labelsize=10)
    plt.tight_layout()
    plt.savefig(path, dpi=300, format='jpg')

    plt.clf()
    plt.close('all')
    gc.collect()

    return



################# RUN #######################

##### produce the plot for one protein  #####

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="MSA representation plot.")
    parser.add_argument("--msa", required=True, help="Path to the MSA.")
    args = parser.parse_args()

    if not args.msa:
        print("Error: The '--msa' argument is missing.")
        parser.print_help()
    else:
        protein_id = args.msa.split('/')[-1][3:-6]

        if os.path.exists(args.msa):
            if args.msa[-5:]!='fasta': print('your file should be in fasta format'); sys.exit(1)
            pathToSave = '/'.join(args.msa.split('/')[:-1]) + f'/{protein_id}_msaRepresentation.jpg'
            #if os.path.isfile(pathToSave):print('done')
            plot_msa_v2(args.msa,d_aa_id, pathToSave, ' ', sort_lines=True, dpi=300)
        else:
            print('The provided fasta file does not exist')




