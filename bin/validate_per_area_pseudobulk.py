#!/usr/bin/env python3
"""
validate_per_area_pseudobulk.py

This version accounts for scale differences:
 - RCTD values are in log scale → exponentiated before aggregation
 - cell2location values are in linear scale (proportion of expression)
All are normalized by total sum per sample–area for scale comparability.

Directory conventions:
 - RCTD:              {rctd_base}/{sample}/gene_expr_ct/{area}/{sample}_{celltype}_gene_expr.csv
 - cell2location_map: {cell2loc_base}/{sample}/gene_expr_ct_mean/{area}/{sample}_{celltype}_gene_expr.csv
"""

import argparse
import os
import sys
import glob
import numpy as np
import pandas as pd
from scipy import stats
import math
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics.pairwise import cosine_similarity

try:
    import anndata
except Exception:
    anndata = None


# ---------- HELPER FUNCTIONS ----------

def read_visium_pseudobulk(h5ad_path, area_label):
    """Compute mean expression per gene for Visium spots in the given area."""
    if not os.path.exists(h5ad_path):
        raise FileNotFoundError(f"h5ad not found: {h5ad_path}")

    if anndata is None:
        raise ImportError("anndata not installed in environment.")

    adata = anndata.read_h5ad(h5ad_path)
    if 'manual_delineation' not in adata.obs:
        raise KeyError(f"obs['manual_delineation'] not found in {h5ad_path}")

    obs_mask = adata.obs['manual_delineation'].astype(str) == str(area_label)
    if obs_mask.sum() == 0:
        print(f"⚠️ No Visium spots found for area '{area_label}' in {h5ad_path}", file=sys.stderr)
        return pd.Series(dtype=float)

    X = adata[obs_mask].X
    if hasattr(X, "toarray"):
        X = X.toarray()

    mean_expr = np.asarray(X).mean(axis=0)
    genes = adata.var_names.tolist()

    series = pd.Series(mean_expr, index=genes, name='visium_mean')

    # Normalize to comparable total sum (CPM-like)
    series = series / series.sum() * 1e6
    return series


def read_method_csvs(base_dir, sample, area, method, masked_celltypes=None):
    """
    Read per-celltype gene expression for a method, aggregate across celltypes.

    For RCTD:   read from gene_expr_ct/{area}/
    For cell2loc: read only from gene_expr_ct_mean/{area}/
    """
    if method.lower() == "rctd":
        search_dir = os.path.join(base_dir, sample, "gene_expr_ct", area)
    elif method.lower() == "cell2loc":
        search_dir = os.path.join(base_dir, sample, "gene_expr_ct_mean", area)
    else:
        raise ValueError("method must be 'rctd' or 'cell2loc'")

    if not os.path.isdir(search_dir):
        print(f"⚠️ Directory not found: {search_dir}", file=sys.stderr)
        return {}, pd.Series(dtype=float)

    files = sorted(glob.glob(os.path.join(search_dir, f"{sample}_*_gene_expr.csv")))
    if not files:
        print(f"⚠️ No gene expression CSVs found in {search_dir}", file=sys.stderr)
        return {}, pd.Series(dtype=float)

    celltype_series = {}

    for fpath in files:
        fname = os.path.basename(fpath)
        celltype = fname.replace(f"{sample}_", "").replace("_gene_expr.csv", "")

        if masked_celltypes:
            masked = [m.strip() for m in masked_celltypes.split(",") if m.strip()]
            if celltype in masked:
                continue

        try:
            df = pd.read_csv(fpath, index_col=0, header=0)
            if df.shape[1] == 1:
                series = df.iloc[:, 0].astype(float)
            else:
                numeric_cols = df.select_dtypes(include=[np.number]).columns
                series = df[numeric_cols[0]].astype(float)
        except Exception as e:
            print(f"⚠️ Could not read {fpath}: {e}", file=sys.stderr)
            continue

        # Convert RCTD log scale to linear
        if method.lower() == "rctd":
            series = np.exp(series)

        series.index = series.index.astype(str)
        celltype_series[celltype] = series

    if not celltype_series:
        return {}, pd.Series(dtype=float)

    # Aggregate across celltypes (sum)
    df_all = pd.concat(celltype_series.values(), axis=1).fillna(0.0)
    agg = df_all.sum(axis=1)

    # Normalize to total sum for scale comparability
    if agg.sum() > 0:
        agg = agg / agg.sum() * 1e6

    agg.name = f"{method}_aggregated"
    return celltype_series, agg


def compute_metrics(true_s, pred_s):
    df = pd.concat([true_s, pred_s], axis=1, join='inner').dropna()
    if df.empty:
        return dict(pearson_r=np.nan, spearman_r=np.nan, mae=np.nan, rmse=np.nan, n_genes=0)

    a, b = df.iloc[:, 0], df.iloc[:, 1]
    pearson_r, _ = stats.pearsonr(a, b) if len(df) > 2 else (np.nan, np.nan)
    spearman_r, _ = stats.spearmanr(a, b) if len(df) > 2 else (np.nan, np.nan)
    mae = np.mean(np.abs(a - b))
    rmse = math.sqrt(np.mean((a - b) ** 2))
    return dict(pearson_r=pearson_r, spearman_r=spearman_r, mae=mae, rmse=rmse, n_genes=len(df))


def plot_scatter(true_s, pred_s, out_png, title=None):
    df = pd.concat([true_s, pred_s], axis=1, join='inner').dropna()
    if df.empty:
        return
    plt.figure(figsize=(5, 5))
    plt.scatter(df.iloc[:, 0], df.iloc[:, 1], s=5, alpha=0.6)
    lims = [df.min().min(), df.max().max()]
    plt.plot(lims, lims, 'k--', lw=1)
    plt.xlabel("Visium (normalized)")
    plt.ylabel("Estimated (normalized)")
    if title:
        plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()


# ---------- MAIN ----------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--area", required=True)
    ap.add_argument("--samples", required=True, help="Comma-separated list of samples")
    ap.add_argument("--spatial-input", required=True)
    ap.add_argument("--rctd-base", required=True)
    ap.add_argument("--cell2loc-base", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--masked-celltypes", default="")
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    samples = [s.strip() for s in args.samples.split(",") if s.strip()]

    metrics_records = []
    scatter_data = []
    summary = pd.DataFrame()   # <-- FIX: ensure summary always exists

    for sample in samples:
        print(f"=== Processing {sample} / {args.area} ===")
        h5ad_path = os.path.join(args.spatial_input, sample, "outs", "matrices",
                                 f"adata_{sample}_manual_delineation_good.h5ad")
        if not os.path.exists(h5ad_path):
            h5ad_path = os.path.join(
                args.spatial_input, sample, "outs", "matrices",
                f"{sample}_manual_delineation.h5ad"
            )
        try:
            visium_series = read_visium_pseudobulk(h5ad_path, args.area)
        except Exception as e:
            print(f"⚠️ Visium read failed for {sample}: {e}", file=sys.stderr)
            visium_series = pd.Series(dtype=float)

        # Read both methods
        _, rctd_agg = read_method_csvs(args.rctd_base, sample, args.area, "rctd", args.masked_celltypes)
        _, c2l_agg = read_method_csvs(args.cell2loc_base, sample, args.area, "cell2loc", args.masked_celltypes)

        sample_out = os.path.join(args.output_dir, sample)
        os.makedirs(sample_out, exist_ok=True)

        if not visium_series.empty:
            visium_series.to_csv(os.path.join(sample_out, f"{sample}_{args.area}_visium_norm.csv"))
        if not rctd_agg.empty:
            rctd_agg.to_csv(os.path.join(sample_out, f"{sample}_{args.area}_rctd_norm.csv"))
        if not c2l_agg.empty:
            c2l_agg.to_csv(os.path.join(sample_out, f"{sample}_{args.area}_cell2loc_norm.csv"))

        for method, est in [("RCTD", rctd_agg), ("cell2loc", c2l_agg)]:
            if est.empty or visium_series.empty:
                continue
            common_genes = visium_series.index.intersection(est.index)
            if len(common_genes) == 0:
                print(f"⚠️ No common genes for {sample} - {method} - {args.area}; skipping.")
                continue

            true_expr = visium_series.loc[common_genes]
            pred_expr = est.loc[common_genes]
        
            if method == "RCTD":
                pred_expr = np.exp(np.clip(pred_expr, -20, 20))
            # Normalize to comparable scales (optional but helps)
            true_expr = true_expr / true_expr.sum()
            pred_expr = pred_expr / pred_expr.sum()

            m = compute_metrics(visium_series, est)
            scatter_data.append((true_expr.values, pred_expr.values, sample, args.area, method))
        
            m.update(sample=sample, area=args.area, method=method)
            metrics_records.append(m)

            out_png = os.path.join(sample_out, f"{sample}_{args.area}_{method}_scatter.png")
            plot_scatter(visium_series, est, out_png, title=f"{sample} {args.area} {method}")

    if metrics_records:
        df = pd.DataFrame(metrics_records)
        df.to_csv(os.path.join(args.output_dir, f"{args.area}_metrics_per_sample.csv"), index=False)
        summary = df.groupby("method").mean(numeric_only=True).reset_index()
        summary.to_csv(os.path.join(args.output_dir, f"{args.area}_metrics_summary.csv"), index=False)

    print(f"✅ Validation for {args.area} completed: results in {args.output_dir}")

    # === Visualization section ===
    import seaborn as sns
    import matplotlib.pyplot as plt
    from sklearn.metrics.pairwise import cosine_similarity

    fig_dir = os.path.join(args.output_dir, "Plots")
    os.makedirs(fig_dir, exist_ok=True)
    sns.set(style="whitegrid", font_scale=1.2)

    # --- 1. Summary barplots ---
    if not summary.empty:
        fig, axes = plt.subplots(1, 2, figsize=(10, 5))
        sns.barplot(data=summary, x="method", y="pearson_r", hue="method", ax=axes[0], palette="viridis", legend=False)
        sns.barplot(data=summary, x="method", y="spearman_r", hue="method", ax=axes[1], palette="magma", legend=False)
        axes[0].set_title("Pearson correlation")
        axes[1].set_title("Spearman correlation")
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f"{args.area}_summary_correlations.png"))
        plt.close()
    else:
        print(f"⚠️ Skipping plots: summary empty for {args.area}")

    # --- 2. Scatterplots ---
    for entry in scatter_data:
        true_expr, pred_expr, sample, label, method = entry
        plt.figure(figsize=(6,6))
        plt.scatter(true_expr, pred_expr, s=6, alpha=0.4)
        lim = max(true_expr.max(), pred_expr.max())
        plt.plot([0, lim], [0, lim], "r--")
        plt.xlabel("True mean expression")
        plt.ylabel("Predicted mean expression")
        plt.title(f"{method} - {sample} - {label}")
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f"{sample}_{label}_{method}_scatter.png"), dpi=300)
        plt.close()


if __name__ == "__main__":
    main()

