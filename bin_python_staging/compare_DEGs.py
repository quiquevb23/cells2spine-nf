#!/usr/bin/env python3
"""
compare_DEGs.py
Compare differential-expression (DEG) results between Cell2location and CSIDE (RCTD)
for each area and cell type.

Outputs:
  • logFC scatter plots with Spearman correlations
  • Venn diagrams for top-N up/down regulated genes
  • Summary CSV of correlations + overlaps
"""
print("beginning script comparing DEGs")

import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
from matplotlib_venn import venn2
import re

# ---------------------------------------------------------------------
def read_deg_table(path):
    """Read DEG CSV and standardize key columns (gene, logFC, Pval)."""
    df = pd.read_csv(path)
    df.columns = [c.strip().strip('"') for c in df.columns]
    gene_col = [c for c in df.columns if "gene" in c.lower()][0]
    logfc_col = [c for c in df.columns if "logfc" in c.lower()][0]
    p_cols = [c for c in df.columns if "p.val" in c.lower() or "p.value" in c.lower()]
    p_col = p_cols[0] if p_cols else None
    df = df[[gene_col, logfc_col] + ([p_col] if p_col else [])].dropna()
    df.columns = ["gene", "logFC"] + (["Pval"] if p_col else [])
    #if "Pval" in df.columns:
    #    df["rank"] = df["Pval"].rank(method="dense", ascending=True)
    #else:
    # fallback: rank by |logFC|
    df["rank"] = df["logFC"].rank(method="dense", ascending=False)
    return df

# ---------------------------------------------------------------------
def top_genes_by_logfc(df, n=100):
    """Return top n up and down genes by logFC."""
    up = df.sort_values("logFC", ascending=False).head(n)["gene"]
    down = df.sort_values("logFC", ascending=True).head(n)["gene"]
    return set(up), set(down)

# ---------------------------------------------------------------------
def compare_methods(df1, df2, area, ct, outdir, n_top=100):
    """Compare two DEG tables: correlations + overlaps + plots."""
    merged = pd.merge(df1, df2, on="gene", suffixes=("_c2l", "_cside"))
    if merged.empty:
        print(f"[WARN] No overlapping genes for {area}/{ct}")
        return None

    # --- correlations ---
    rho_rank, _ = spearmanr(merged["rank_c2l"], merged["rank_cside"])
    rho_fc, _ = spearmanr(merged["logFC_c2l"], merged["logFC_cside"])

    # --- scatter plot ---
    plt.figure(figsize=(5,5))
    sns.scatterplot(
        x="logFC_c2l", y="logFC_cside", data=merged, s=18, alpha=0.6, edgecolor=None
    )
    plt.axhline(0, color="gray", lw=0.5)
    plt.axvline(0, color="gray", lw=0.5)
    plt.title(f"{area} – {ct}\nSpearman ρ(logFC)={rho_fc:.2f}, ρ(rank)={rho_rank:.2f}")
    plt.xlabel("Cell2location logFC")
    plt.ylabel("CSIDE logFC")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"{ct}_logFC_scatter.png"), dpi=300)
    plt.close()

    # --- top-N up/down overlaps ---
    up1, down1 = top_genes_by_logfc(df1, n_top)
    up2, down2 = top_genes_by_logfc(df2, n_top)

    overlaps = {}
    for label, s1, s2 in [("up", up1, up2), ("down", down1, down2)]:
        overlap = len(s1 & s2)
        overlaps[f"{label}_overlap"] = overlap
        plt.figure(figsize=(4,4))
        venn2([s1, s2], set_labels=("Cell2location", "CSIDE"))
        plt.title(f"{area} – {ct} ({label}regulated)\nTop {n_top} overlap={overlap}")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"{ct}_{label}_top{n_top}_venn.png"), dpi=300)
        plt.close()

    return {
        "area": area,
        "celltype": ct,
        "rho_rank": rho_rank,
        "rho_logFC": rho_fc,
        "up_overlap": overlaps["up_overlap"],
        "down_overlap": overlaps["down_overlap"],
        "n_genes_shared": len(merged)
    }

# ---------------------------------------------------------------------
def main(args):
    area = args.area
    cell2loc_base = args.cell2loc_dir
    rctd_base = args.rctd_dir
    out_base = args.output_dir
    os.makedirs(out_base, exist_ok=True)
    top_n = args.top_n

    c2l_dir = os.path.join(cell2loc_base, area)
    rctd_dir = os.path.join(rctd_base, area)

    if not os.path.isdir(rctd_dir):
        print(f"[WARN] Missing CSIDE DEGs for {area}")
        return
    if not os.path.isdir(c2l_dir):
        print(f"[WARN] Missing Cell2location DEGs for {area}")
        return
    print(rctd_dir)
    # detect cell types by matching filenames
    cts = []
    for f in os.listdir(rctd_dir):
        if f.startswith("DEA_") and f.endswith("_limma.csv"):
            m = re.match(r"DEA_(.+)_(.+)_limma\.csv", f)
            if m:
                area_in_file, ct = m.groups()
                if area_in_file == area:
                    cts.append(ct)
    cts = [
        ct for ct in cts
        if os.path.exists(os.path.join(c2l_dir, f"DEA_{area}_{ct}_cell2loc_limma.csv"))
    ]
    print(f"cell types detected for area {area}:\n {cts}")
    summary_rows = []
    for ct in cts:
        f1 = os.path.join(c2l_dir, f"DEA_{area}_{ct}_cell2loc_limma.csv")
        f2 = os.path.join(rctd_dir, f"DEA_{area}_{ct}_limma.csv")
        df1, df2 = read_deg_table(f1), read_deg_table(f2)
        res = compare_methods(df1, df2, area, ct, out_base, n_top=top_n)
        if res:
            summary_rows.append(res)

    if summary_rows:
        pd.DataFrame(summary_rows).to_csv(
            os.path.join(out_base, f"DEG_overlap_summary_{area}.csv"), index=False
        )
        print(f"[INFO] Wrote summary for {len(summary_rows)} cell types -> {out_base}")

# ---------------------------------------------------------------------
if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Compare DEGs between Cell2location and CSIDE")
    p.add_argument("--area", required=True)
    p.add_argument("--cell2loc_dir", required=True)
    p.add_argument("--rctd_dir", required=True)
    p.add_argument("--output_dir", required=True)
    p.add_argument("--top_n", type=int, default=1000, help="Top N genes for overlap plots")
    args = p.parse_args()
    main(args)
