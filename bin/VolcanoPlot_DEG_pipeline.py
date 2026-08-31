"""
Generate volcano plots for DEA areas
"""

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from adjustText import adjust_text


def parse_args():
    parser = argparse.ArgumentParser(description="Volcano plots for DEA areas")
    parser.add_argument("--output_base_dir", required=True)
    parser.add_argument("--samples", nargs="+", required=True)
    return parser.parse_args()


args = parse_args()

output_dir = os.path.join(args.output_base_dir, "DEA_areas_good")
os.makedirs(output_dir, exist_ok=True)

dirs = [
    d for d in os.listdir(output_dir)
    if os.path.isdir(os.path.join(output_dir, d))
    and not d.startswith("Spatial")
]

for area in dirs:
    area_dir = os.path.join(output_dir, area)
    plot_dir = os.path.join(area_dir, "ViolinPlots_DEGs")
    os.makedirs(plot_dir, exist_ok=True)

    for f in os.listdir(area_dir):
        if not f.endswith("pseudobulk.csv"):
            continue

        df = pd.read_csv(os.path.join(area_dir, f))
        
        if "adj.P.Val" in df.columns and "FDR" not in df.columns:
            df = df.rename(columns={"adj.P.Val": "FDR"})

        degs = df[df["FDR"] < 0.05]

        up = degs[degs["logFC"] > 0]
        down = degs[degs["logFC"] < 0]

        fig, ax = plt.subplots(figsize=(7, 5))

        ax.scatter(df["logFC"], -np.log10(df["FDR"]), s=1, alpha=0.4)
        ax.scatter(up["logFC"], -np.log10(up["FDR"]), s=3, color="red")
        ax.scatter(down["logFC"], -np.log10(down["FDR"]), s=3, color="blue")

        ax.axvline(0, linestyle="--")
        ax.axhline(1.301, linestyle="--")

        texts = []
        for _, r in up.nlargest(10, "logFC").iterrows():
            texts.append(ax.text(r["logFC"], -np.log10(r["FDR"]), r["gene"], fontsize=8))
        for _, r in down.nsmallest(10, "logFC").iterrows():
            texts.append(ax.text(r["logFC"], -np.log10(r["FDR"]), r["gene"], fontsize=8))

        adjust_text(texts)

        ax.set_xlabel("log₂FC")
        ax.set_ylabel("-log₁₀(FDR)")
        ax.set_title(f"Volcano plot – {area}")

        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, f"{area}_volcano_plot.png"), dpi=300)
        plt.close()

