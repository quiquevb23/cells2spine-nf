#!/usr/bin/env python3
"""
plot_global_FE_heatmaps.py

Aggregate the per-area Reactome FE overlap summaries (from compare_FE.py)
and plot two global heatmaps (UP and DOWN Jaccard agreement) for GSEA results.

Usage:
  python3 plot_global_FE_heatmaps.py \
    --input_base <base_dir_with_all_area_subdirs> \
    --output_dir <output_dir_for_heatmaps>
"""

import argparse
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def main(args):
    input_base = args.input_base
    out_dir = args.output_dir
    os.makedirs(out_dir, exist_ok=True)

    # find all per-area summary CSVs recursively
    csv_files = glob.glob(os.path.join(input_base, "*", "FE_GO_BP_overlap_summary_*.csv"))
    if not csv_files:
        print(f"[ERROR] No FE summary CSVs found under {input_base}")
        return

    dfs = []
    for f in csv_files:
        try:
            df = pd.read_csv(f)
            dfs.append(df)
        except Exception as e:
            print(f"[WARN] Failed to read {f}: {e}")

    if not dfs:
        print("[ERROR] No valid summary data loaded.")
        return

    all_df = pd.concat(dfs, ignore_index=True)
    all_df = all_df.query('analysis == "GSEA"').copy()

    if all_df.empty:
        print("[WARN] No GSEA entries found.")
        return

    # --- Build heatmap matrices ---
    up_pivot = all_df.pivot(index="area", columns="celltype", values="up_jaccard")
    down_pivot = all_df.pivot(index="area", columns="celltype", values="down_jaccard")

    # sort by area and celltype (optional)
    up_pivot = up_pivot.sort_index(axis=0)
    up_pivot = up_pivot.reindex(sorted(up_pivot.columns), axis=1)
    down_pivot = down_pivot.loc[up_pivot.index, up_pivot.columns]

    # --- Plot helper ---
    def plot_heatmap(data, title, fname):
        plt.figure(figsize=(1.2 * len(data.columns) + 3, 0.8 * len(data.index) + 3))
        sns.heatmap(
            data,
            annot=True,
            fmt=".2f",
            cmap="viridis",
            linewidths=0.5,
            linecolor="gray",
            cbar_kws={"label": "Jaccard overlap"},
            vmin=0,
            vmax=1,
            mask=data.isna()
        )
        plt.title(title, fontsize=14, pad=12)
        plt.xlabel("Cell type")
        plt.ylabel("Area")
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, fname), dpi=300)
        plt.close()

    # --- Plot both ---
    plot_heatmap(up_pivot, "GSEA GO_BP UP overlap (Cell2location vs CSIDE)", "FE_GSEA_GO_BP_UP_heatmap.png")
    plot_heatmap(down_pivot, "GSEA GO_BP DOWN overlap (Cell2location vs CSIDE)", "FE_GSEA_GO_BP_DOWN_heatmap.png")

    print(f"[INFO] Saved heatmaps to {out_dir}")

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Plot global FE GO_BP agreement heatmaps for GSEA")
    p.add_argument("--input_base", required=True, help="Base directory where area subfolders contain FE_GO_BP_overlap_summary_<area>.csv")
    p.add_argument("--output_dir", required=True, help="Directory to save the heatmaps")
    args = p.parse_args()
    main(args)
