import pandas as pd
import seaborn as sns
import numpy as np
from pathlib import Path
import os
import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import FormatStrFormatter


WORK_DIR = Path.cwd()

#%% Segmentation plot
def plot_segmentation(dir_in, segments, ptm=pd.DataFrame(), width=20, height=5):

    segments_grouped = segments.groupby('protein')

    for protein, group in segments_grouped:

        data = pd.read_csv(os.path.join(dir_in, f"{protein}_GEMME_pLDDT.csv"))
        col = 'GEMME_mean'
        mu_mat = data[col].values
        n_res = len(mu_mat)

        signal = pd.DataFrame({
            'x': list(range(1, len(data) + 1)),
            'y': list(data[col]),
            'type': ['GEMME'] * len(data)
        })

        gemme = group[group['type'] == 'GEMME']
        cp = pd.DataFrame({
            'x': list(gemme['end'][:-1]),'type': ['GEMME'] * (len(gemme) - 1)})
        mean = pd.DataFrame({
            'x': list(gemme['start'] - 1),'xend': list(gemme['end']),
            'y': list(gemme['mean']),'type': ['GEMME'] * len(gemme),
            'state': list(gemme['state'])
        })
        
        if 'pLDDT' in data.columns:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(width, height), gridspec_kw={'height_ratios': [1, 6]})
            # Plot GEMME on the top subplot
            for _, row in mean[mean['type'] == 'GEMME'].iterrows():
                ax2.axvspan(row['x'], row['xend'], color='white' if row['state'] == 0 else 'red' if row['state'] == 1 else 'purple', alpha=0.1)
            ax2.plot(signal[signal['type'] == 'GEMME']['x'], signal[signal['type'] == 'GEMME']['y'], label='GEMME', linewidth=1.3)
            for _, row in cp[cp['type'] == 'GEMME'].iterrows():
                ax2.axvline(row['x'], color='red', alpha=0.9, linewidth=1.3)
            for _, row in mean[mean['type'] == 'GEMME'].iterrows():
                ax2.plot([row['x'], row['xend']], [row['y'], row['y']], color='black',linewidth=1, alpha=0.75)
            ax2.set_ylabel('Mutational sensitivity', fontsize=19, labelpad=10)
            

            # Plot pLDDT heatmap on the bottom subplot
            palette_bis={0:'#FF7D45', 1:'#FFDB13', 2:'#65CBF3', 3:'#0053D6'}
            bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
            norm = BoundaryNorm(bounds, len(palette_bis))
            bins = [0, 50, 70, 90, 100]
            df = pd.DataFrame(list(data['pLDDT']*100), columns=['value'])
            # Bin the values into groups
            df['group'] = pd.cut(df['value'], bins=bins, labels=[0,1,2,3], right=False)
            grouped_values_numeric = np.array(df['group'])

            # Reshape to a one-row matrix with 395 columns
            data_matrix = grouped_values_numeric.reshape(1, -1)

            # Map the palette to colors for the heatmap
            colors = [palette_bis[val] for val in sorted(palette_bis.keys())]
            cmap = sns.color_palette(colors, len(colors))
            sns.heatmap(data_matrix, cmap=cmap,norm=norm, linewidths=0.2, cbar=False, xticklabels=False, yticklabels=False, alpha=0.85, ax=ax1)

            # Adjust xticks for ax2
            xticks = range(50, len(data) + 1, 50)
            
            ax2.set_xticks(xticks)
            ax2.set_xticklabels([str(x) for x in xticks], fontsize=15)
            ax2.set_ylim(0, 1)
            #ax2.set_yticklabels([str(round(x, 1)) for x in ax2.get_yticks()], fontsize=15)
            ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            ax2.tick_params(axis='y', labelsize=15)
    
            # Remove x-axis ticks and labels from ax1
            ax1.set_xticks([])
            ax1.set_xticklabels([])
            ax1.set_xlim(0, len(data))
            ax2.set_xlim(1, len(data))

            
            plt.xlabel('Residue index', fontsize=20, labelpad=12)
            plt.subplots_adjust(hspace=0.05)
            plt.tight_layout()
            output = f'{dir_in}/{protein}_SegProfile.png'
            plt.savefig(output, dpi=700)
            plt.close(fig)



        else:
            fig,ax2 = plt.subplots(1, 1, figsize=(width, height))
            # Plot GEMME on the top subplot
            for _, row in mean[mean['type'] == 'GEMME'].iterrows():
                ax2.axvspan(row['x'], row['xend'], color='white' if row['state'] == 0 else 'red' if row['state'] == 1 else 'purple', alpha=0.1)
            ax2.plot(signal[signal['type'] == 'GEMME']['x'], signal[signal['type'] == 'GEMME']['y'], label='GEMME', linewidth=1.3)
            for _, row in cp[cp['type'] == 'GEMME'].iterrows():
                ax2.axvline(row['x'], color='red', alpha=0.9, linewidth=1.3)
            for _, row in mean[mean['type'] == 'GEMME'].iterrows():
                ax2.plot([row['x'], row['xend']], [row['y'], row['y']], color='black',linewidth=1, alpha=0.75)
            ax2.set_ylabel('Mutational sensitivity', fontsize=19, labelpad=10)
    

            # Adjust xticks for ax2
            xticks = range(50, len(data) + 1, 50)
            
            ax2.set_xticks(xticks)
            ax2.set_xticklabels([str(x) for x in xticks], fontsize=15)
            ax2.set_ylim(0, 1)
            ax2.set_xlim(1, len(data))
            ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            ax2.tick_params(axis='y', labelsize=15) 
            
        
            plt.xlabel('Residue index', fontsize=20, labelpad=12)
            plt.tight_layout()
            output = dir_in+f'/{protein}_SegProfile.png'
            plt.savefig(output, dpi=700)
            plt.close(fig)


if __name__ == '__main__':
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Segmentation plot.")
    
    # Long-style arguments with -- prefix
    parser.add_argument("--seg-file", required=True, help="Path to the input segmentation file.")
    parser.add_argument("--input-dir", default='', help="Path to the file with PTMs data to plot.")
    parser.add_argument("--ptm-file", default='', help="Path to the file with PTMs data to plot.")

    # Parse arguments
    args = parser.parse_args()


    # Check if all required arguments are provided
    if not args.seg_file:
        print("Error: The '--seg-file' argument is missing.")
        parser.print_help()
    elif not args.input_dir:
        print("Error: The '--input-dir' argument is missing.")
        parser.print_help()
    else:
        # Run the function with the command-line arguments
        
        print(f'Generating segmentation plots from {args.seg_file} data')

        segments = pd.read_csv(f'{args.seg_file}')

        if args.ptm_file:
            print('PTM file is provided')
            df_ptm = pd.read_csv(f'{WORK_DIR}/{args.ptm_file}')
            plot_segmentation(args.input_dir, segments, ptm=df_ptm)
        else:
            plot_segmentation(args.input_dir, segments)
