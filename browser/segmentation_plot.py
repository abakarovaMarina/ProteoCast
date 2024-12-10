import pandas as pd
import numpy as np
from pathlib import Path
from utils import read_GEMME_mat, bfactorsPDB
import os
from Bio import SeqIO
import gemmi as gm
import argparse
import matplotlib.pyplot as plt

WORK_DIR = Path.cwd()

#%% Segmentation plot
def plot_segmentation(dir_prot, dir_save, segments, ptm=pd.DataFrame(), width=20, height=8):

    segments_grouped = segments.groupby('protein')

    for protein, group in segments_grouped:
        data = pd.read_csv(os.path.join(WORK_DIR, dir_prot, f"{protein}_sensitive_residues.csv"))
        print('PROTEIN:', protein)
        signal = pd.DataFrame({
            'x': list(range(1, len(data) + 1)) * (2 if 'pLDDT' in data else 1),
            'y': list(data['GEMME_mean']) + (list(data['pLDDT']) if 'pLDDT' in data else []),
            'type': ['GEMME'] * len(data) + (['pLDDT'] * len(data) if 'pLDDT' in data else [])
        })

        gemme = group[group['type'] == 'GEMME']
        plddt = group[group['type'] == 'pLDDT'] if 'pLDDT' in data else pd.DataFrame()
        
        cp = pd.DataFrame({
            'x': list(gemme['end'][:-1]) + (list(plddt['end'][:-1]) if not plddt.empty else []),
            'type': ['GEMME'] * (len(gemme) - 1) + (['pLDDT'] * (len(plddt) - 1) if not plddt.empty else [])
        })
        
        mean = pd.DataFrame({
            'x': list(gemme['start'] - 1) + (list(plddt['start'] - 1) if not plddt.empty else []),
            'xend': list(gemme['end']) + (list(plddt['end']) if not plddt.empty else []),
            'y': list(gemme['mean']) + (list(plddt['mean']) if not plddt.empty else []),
            'type': ['GEMME'] * len(gemme) + (['pLDDT'] * len(plddt) if not plddt.empty else []),
            'state': list(gemme['state']) + (list(plddt['state']) if not plddt.empty else [])
        })
        
        if 'pLDDT' not in data.columns or plddt.empty:
            # Single plot case
            height = 4
            fig, ax = plt.subplots(1, 1, figsize=(width, height))
            
            for _, row in mean[mean['type'] == 'GEMME'].iterrows():
                ax.axvspan(row['x'], row['xend'], color='gold' if row['state'] == 0 else 'red' if row['state'] == 1 else 'purple', alpha=0.1)
            ax.plot(signal[signal['type'] == 'GEMME']['x'], signal[signal['type'] == 'GEMME']['y'], label='GEMME')
            for _, row in cp[cp['type'] == 'GEMME'].iterrows():
                ax.axvline(row['x'], color='red')
            for _, row in mean[mean['type'] == 'GEMME'].iterrows():
                ax.plot([row['x'], row['xend']], [row['y'], row['y']], color='blue')
            
            ax.set_ylabel('Average GEMME score', fontsize=17, labelpad=10)
            ax.set_xlim(1, len(data))
            xticks = range(100, len(data) + 1, 100) if len(data) > 1500 else range(50, len(data) + 1, 50)
            ax.set_xticks(xticks)
            ax.set_xticklabels([str(x) for x in xticks], fontsize=15)
            ax.yaxis.set_major_locator(plt.FixedLocator(ax.get_yticks()))
            ax.set_yticklabels([str(round(x, 1)) for x in ax.get_yticks()], fontsize=15)


            for _, row in ptm.iterrows():
                ax.axvline(int(row['position']), linestyle='--', color='lightgrey')
            
            plt.xlabel('Residue index', fontsize=20, labelpad=12)
            plt.tight_layout()
            output = f'{WORK_DIR}/{dir_save}/{protein}_segmentation_gemme_only.png'
            plt.savefig(output, dpi=500)
            plt.close(fig)
            
    
        else:
            # Dual plot case
            height = 8
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(width, height))
            
            # GEMME Plot
            for _, row in mean[mean['type'] == 'GEMME'].iterrows():
                ax2.axvspan(row['x'], row['xend'], color='gold' if row['state'] == 0 else 'red' if row['state'] == 1 else 'purple', alpha=0.1)
            ax2.plot(signal[signal['type'] == 'GEMME']['x'], signal[signal['type'] == 'GEMME']['y'], label='GEMME')
            for _, row in cp[cp['type'] == 'GEMME'].iterrows():
                ax2.axvline(row['x'], color='red')
            for _, row in mean[mean['type'] == 'GEMME'].iterrows():
                ax2.plot([row['x'], row['xend']], [row['y'], row['y']], color='blue')
            
            ax2.set_ylabel('Average GEMME score', fontsize=17, labelpad=10)
            
            # pLDDT Plot
            for _, row in mean[mean['type'] == 'pLDDT'].iterrows():
                ax1.axvspan(row['x'], row['xend'], color='gold' if row['state'] == 0 else 'red' if row['state'] == 1 else 'purple', alpha=0.1)
            ax1.plot(signal[signal['type'] == 'pLDDT']['x'], signal[signal['type'] == 'pLDDT']['y'], label='pLDDT')
            for _, row in cp[cp['type'] == 'pLDDT'].iterrows():
                ax1.axvline(row['x'], color='red')
            for _, row in mean[mean['type'] == 'pLDDT'].iterrows():
                ax1.plot([row['x'], row['xend']], [row['y'], row['y']], color='blue')
            
            ax1.set_ylabel('pLDDT score', fontsize=17, labelpad=10)
            
            # Configure x-axis
            xticks = range(100, len(data) + 1, 100) if len(data) > 1500 else range(50, len(data) + 1, 50)
            ax2.set_xlim(1, len(data))
            ax1.set_xlim(1, len(data))
            ax2.set_xticks(xticks)
            ax2.set_xticklabels([str(x) for x in xticks], fontsize=15)
            ax1.set_xticks([])
            ax1.set_xticklabels([])
            
            # y-axis formatting
            ax1.yaxis.set_major_locator(plt.FixedLocator(ax1.get_yticks()))
            ax1.set_yticklabels([str(round(x, 1)) for x in ax1.get_yticks()], fontsize=15)
            ax2.yaxis.set_major_locator(plt.FixedLocator(ax2.get_yticks()))
            ax2.set_yticklabels([str(round(x, 1)) for x in ax2.get_yticks()], fontsize=15)
            
            for _, row in ptm.iterrows():
                ax2.axvline(int(row['position']), linestyle='--', color='lightgrey')
            
            plt.xlabel('Residue index', fontsize=20, labelpad=12)
            plt.subplots_adjust(hspace=0.05)
            plt.tight_layout()
            output = f'{WORK_DIR}/{dir_save}/{protein}_segmentation.png'
            plt.savefig(output, dpi=500)
            plt.close(fig)

if __name__ == '__main__':
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Split a FASTA file into smaller chunks.")
    
    # Long-style arguments with -- prefix
    parser.add_argument("--input", required=True, help="Path to the input segmentation file.")
    parser.add_argument("--input-dir", default="csv", help="Directory to the input GEMME and plDDT signals.")
    parser.add_argument("--output-dir", default="segmentation", help="Directory to save the split FASTA files.")
    parser.add_argument("--ptm-file", default='', help="Path to the file with PTMs data to plot.")

    # Parse arguments
    args = parser.parse_args()


    # Check if all required arguments are provided
    if not args.input:
        print("Error: The '--input' argument is missing.")
        parser.print_help()
    elif not args.input_dir:
        print("Error: The '--input-dir' argument is missing.")
        parser.print_help()
    elif not args.output_dir:
        print("Error: The '--output-dir' argument is missing.")
        parser.print_help()
    else:
        # Run the function with the command-line arguments
        
        print(f'Generating segmentation plots from {args.input} data')
        print(f'Output directory: {args.output_dir}')
        print(f'Input directory: {args.input_dir}')
        print(f'input file : {args.input}')

        segments = pd.read_csv(f'{WORK_DIR}/{args.input}')

        if not os.path.exists(f'{WORK_DIR}/{args.output_dir}/'):
            print(f"Warning: The directory {args.output_dir} does not exist. Creating it.")
            os.makedirs(f'{WORK_DIR}/{args.output_dir}/')

        if args.ptm_file:
            print('PTM file is provided')
            df_ptm = pd.read_csv(f'{WORK_DIR}/{args.ptm_file}')
            plot_segmentation(args.input_dir, args.output_dir, segments, ptm=df_ptm)
        else:
            plot_segmentation(args.input_dir, args.output_dir, segments)
