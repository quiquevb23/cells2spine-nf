# Script to generate volcano plots for DEGs 

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
from adjustText import adjust_text

def parse_args():
    parser = argparse.ArgumentParser(description="Process directories for single-cell data.")
    
    # Define the arguments for base_dir and output_base_dir
    parser.add_argument('--base_dir', type=str, required=True, 
                        help="Base directory for input data.")
    parser.add_argument('--output_base_dir', type=str, required=True, 
                        help="Base directory for output data.")
    
    # Parse the arguments
    return parser.parse_args()

args = parse_args()

#base_dir = os.path.join(args.base_dir, "DEA_dorsal")
output_dir = os.path.join(args.output_base_dir, "DEA_areas_good")

#os.makedirs(base_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

# Filter directories inside output_dir excluding those starting with 'Spatial'
dirs = [d for d in os.listdir(output_dir)
        if os.path.isdir(os.path.join(output_dir, d)) and not d.startswith("Spatial")]

for dir in dirs:
    pathdir = os.path.join(output_dir, dir)
    
    # Make output subdirectory for plots
    plot_dir = os.path.join(pathdir, "ViolinPlots_DEGs")
    os.makedirs(plot_dir, exist_ok=True)

    # Process only CSVs ending with 'pseudobulk.csv'
    for f in os.listdir(pathdir):
        if f.endswith("pseudobulk.csv"):
            df = pd.read_csv(os.path.join(pathdir, f))

            # Define differentially expressed genes (DEGs)
            degs = df[df['FDR'] < 0.05]
            nre = len(degs)

            # Upregulated (log2FC > 0 & FDR < 0.05)
            up = degs[degs['logFC'] > 0]
            nre_up = len(up)
            up_top10 = up.nlargest(10, 'logFC')

            # Downregulated (log2FC < 0 & FDR < 0.05)
            down = degs[degs['logFC'] < 0]
            nre_down = len(down)
            down_top10 = down.nsmallest(10, 'logFC')

            min_up_logFC = round(up['logFC'].min(), 3) if not up.empty else None
            max_down_logFC = round(down['logFC'].max(), 3) if not down.empty else None

            # Create figure and axis
            fig, ax = plt.subplots(figsize=(7, 5))

            # Volcano plot
            ax.scatter(df['logFC'], -np.log10(df['FDR']), s=1, color="gray", alpha=0.5)
            ax.scatter(up['logFC'], -np.log10(up['FDR']), s=3, color="red", label="Upregulated")
            ax.scatter(down['logFC'], -np.log10(down['FDR']), s=3, color="blue", label="Downregulated")

            # Cutoff lines
            ax.axvline(0, color="grey", linestyle="--")
            ax.axhline(1.301, color="grey", linestyle="--")  # FDR = 0.05

            # Annotations
            texts = []
            for _, r in up_top10.iterrows():
                texts.append(ax.text(x=r['logFC'], y=-np.log10(r['FDR']), s=r['gene'], fontsize=8, color='red'))
            for _, r in down_top10.iterrows():
                texts.append(ax.text(x=r['logFC'], y=-np.log10(r['FDR']), s=r['gene'], fontsize=8, color='blue'))

            adjust_text(
                texts,
                expand_points=(5, 5),
                expand_text=(5, 5),
                force_points=(5, 5),
                force_text=(4, 4),
                avoid_self=True,
                lim=2000,
                only_move={'points': 'y', 'text': 'xy'},
                arrowprops=dict(arrowstyle="-", color='black', lw=0.5, alpha=0.6,
                                connectionstyle="arc3,rad=0.2")
            )

            # Legend and labels
            legend_text = (f"DEGs with FDR<0.05:\nUpregulated genes: {nre_up}\n"
                           f"Downregulated genes: {nre_down}\nMin logFC values: {min_up_logFC}, {max_down_logFC}")
            ax.legend(title=legend_text, loc="upper left", bbox_to_anchor=(1.05, 1), borderaxespad=0.)
            ax.set_xlabel("log₂FC")
            ax.set_ylabel("-log₁₀(FDR)")
            ax.set_title(f"Volcano Plot for {dir}")

            # Save figure
            plot_filename = f"{dir}_volcano_plot.png"
            plt.tight_layout()
            plt.savefig(os.path.join(plot_dir, plot_filename), dpi=300)
            plt.close()
