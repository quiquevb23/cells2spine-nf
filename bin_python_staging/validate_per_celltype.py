#!/usr/bin/env python3


import argparse
import pandas as pd
import numpy as np
import anndata as ad
from scipy.stats import pearsonr, spearmanr
import os
import re
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import cosine_similarity

def compute_metrics(true_vals, pred_vals):
    mask = (~np.isnan(true_vals)) & (~np.isnan(pred_vals))
    if mask.sum() == 0:
        return {"pearson_r": np.nan, "spearman_r": np.nan, "mae": np.nan, "rmse": np.nan, "n_genes": 0}

    pearson_r = pearsonr(true_vals[mask], pred_vals[mask])[0]
    spearman_r = spearmanr(true_vals[mask], pred_vals[mask])[0]
    mae = np.mean(np.abs(true_vals[mask] - pred_vals[mask]))
    rmse = np.sqrt(np.mean((true_vals[mask] - pred_vals[mask])**2))
    return {"pearson_r": pearson_r, "spearman_r": spearman_r, "mae": mae, "rmse": rmse, "n_genes": mask.sum()}


def parse_conditions_map(conditions_map_list):
    mapping = {}
    for pair in conditions_map_list:
        if ":" not in pair:
            continue
        sample, condition = pair.split(":")
        mapping[sample] = condition
    return mapping

def extract_sample_name(filename):
    """
    Extract sample name like 'Spatial_1' or 'Spatial_12' from filename.
    """
    match = re.search(r"(Spatial_\d+|sp\d+)", filename)
    if match:
        return match.group(1)
    return None


def load_and_aggregate(method_dir, celltype):
    """
    Aggregate expression for a given celltype, separately for each sample (across all areas).
    Returns dict {sample: pd.Series(gene_expr)}
    """
    per_sample = defaultdict(list)
    for root, _, files in os.walk(method_dir):
        for f in files:
            if f.endswith(f"{celltype}_gene_expr.csv"):
                sample = extract_sample_name(f)
                if sample is None:
                    print(f"[WARN] Could not extract sample name from filename: {f}")
                    continue
                df = pd.read_csv(os.path.join(root, f), index_col=0)
                if "RCTD" in root:
                    # RCTD is log scale -> convert to exp
                    df.iloc[:, 0] = np.exp(df.iloc[:, 0])
                per_sample[sample].append(df)

    # Aggregate per sample
    result = {}
    for sample, dfs in per_sample.items():
        if len(dfs) == 0:
            continue
        all_expr = pd.concat(dfs, axis=1)
        result[sample] = all_expr.mean(axis=1)
    return result


def main():
    parser = argparse.ArgumentParser(description="Compare per-celltype expression with single-cell reference by condition.")
    parser.add_argument("--ref_level", required=True, help="Level of cell type annotation (l2 or l3)")
    parser.add_argument("--single_cell_ref", required=True, help="Path to single-cell reference h5ad")
    parser.add_argument("--run_base_dir", required=True)
    parser.add_argument("--celltype", required=True)
    parser.add_argument("--spatial_input", required=True)
    parser.add_argument("--ref_label", required=True)
    parser.add_argument("--conditions_map", nargs="+", required=True, help="List of sample:condition mappings")
    args = parser.parse_args()

    # Directories
    cell2loc_dir = os.path.join(args.run_base_dir, "cell2location_map")
    rctd_dir = os.path.join(args.run_base_dir, "RCTD")

    # Parse mapping
    cond_map = parse_conditions_map(args.conditions_map)

    # Load single-cell reference
    print(f"Loading single-cell reference from {args.single_cell_ref}")
    adata = ad.read_h5ad(args.single_cell_ref)
    obs_key = f"celltypes_{args.ref_level}"
    if obs_key not in adata.obs:
        raise ValueError(f"Cell type level {obs_key} not found in adata.obs!")

    # Check condition field
    if "condition" not in adata.obs.columns:
        raise ValueError("Single-cell reference lacks .obs['condition'] field!")

    # Aggregate single-cell means per (condition, celltype)
    sc_means = {}
    for condition in adata.obs["condition"].unique():
        mask = (adata.obs[obs_key] == args.celltype) & (adata.obs["condition"] == condition)
        if mask.sum() == 0:
            continue
        expr = np.array(adata[mask].X.mean(axis=0)).flatten()
        sc_means[condition] = pd.Series(expr, index=adata.var_names)

    results = []
    detailed_expr = []
    scatter_data = []

    for method, dir_ in [("RCTD", rctd_dir), ("cell2loc", cell2loc_dir)]:
        print(f"Processing {method} for {args.celltype}")
        sample_exprs = load_and_aggregate(dir_, args.celltype)
        for sample, expr in sample_exprs.items():
            if sample not in cond_map:
                print(f"Sample {sample} not found in conditions map; skipping.")
                continue

            condition = cond_map[sample]
            if condition not in sc_means:
                print(f"No single-cell reference found for condition {condition}; skipping {sample}.")
                continue

            ref_expr = sc_means[condition]
            common_genes = ref_expr.index.intersection(expr.index)

            true_expr = ref_expr[common_genes]
            pred_expr = expr[common_genes]

            metrics = compute_metrics(ref_expr[common_genes].values, expr[common_genes].values)
            scatter_data.append((true_expr.values, pred_expr.values, sample, args.celltype, method))
            metrics.update({
                "sample": sample,
                "condition": condition,
                "method": method,
                "celltype": args.celltype,
                "n_genes": len(common_genes)
            })
            results.append(metrics)

            detailed_expr.append(pd.DataFrame({
                "gene": common_genes,
                "single_cell": ref_expr[common_genes].values,
                #method: expr[common_genes].values,
                method: pred_expr.values,
                "sample": sample,
                "condition": condition
            }))

    # Output paths
    outdir = os.path.join(args.run_base_dir, "Validation_Celltype")
    os.makedirs(outdir, exist_ok=True)

    per_sample_df = pd.DataFrame(results)
    if per_sample_df.empty:
        print(
            f"[WARN] No valid comparisons for celltype '{args.celltype}'. "
            "Likely absent in SC reference, spatial outputs, or condition overlap."
        )
    return

    per_sample_path = os.path.join(outdir, f"{args.celltype}_metrics_per_sample.csv")
    per_sample_df.to_csv(per_sample_path, index=False)
    print(f"Saved per-sample metrics to {per_sample_path}")


    print("per_sample_df shape:", per_sample_df.shape)
    print("per_sample_df columns:", per_sample_df.columns.tolist())

    # Summary (mean across samples per method)
    summary = (
        per_sample_df.groupby(["method"])
        .agg({"pearson_r": "mean", "spearman_r": "mean", "mae": "mean", "rmse": "mean", "n_genes": "mean"})
        .reset_index()
    )
    summary_path = os.path.join(outdir, f"{args.celltype}_metrics_summary.csv")
    summary.to_csv(summary_path, index=False)
    print(f"Saved summary metrics to {summary_path}")

    # Save combined detailed expression
    if detailed_expr:
        all_expr = pd.concat(detailed_expr, axis=0)
        all_expr.to_csv(os.path.join(outdir, f"{args.celltype}_expression_comparison.csv"), index=False)

    # === Visualization section ===
    import seaborn as sns
    import matplotlib.pyplot as plt
    from sklearn.metrics.pairwise import cosine_similarity

    fig_dir = os.path.join(outdir, "Plots")
    os.makedirs(fig_dir, exist_ok=True)
    sns.set(style="whitegrid", font_scale=1.2)

    # --- 1. Summary barplots ---
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    sns.barplot(data=summary, x="method", y="pearson_r", ax=axes[0], palette="viridis")
    axes[0].set_title("Mean Pearson correlation")
    sns.barplot(data=summary, x="method", y="spearman_r", ax=axes[1], palette="magma")
    axes[1].set_title("Mean Spearman correlation")
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"{args.celltype}_correlation_summary_barplot.png"), dpi=300)
    plt.close()

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

    try:
        if "pred_matrix" in locals() and "true_matrix" in locals():
            celltypes = list(true_matrix.index)
            sim_matrix = cosine_similarity(pred_matrix, true_matrix)
            sim_df = pd.DataFrame(sim_matrix, index=celltypes, columns=celltypes)
            plt.figure(figsize=(10, 8))
            sns.heatmap(sim_df, cmap="coolwarm", annot=False)
            plt.title("Cosine similarity between predicted and true cell-type profiles")
            plt.tight_layout()
            plt.savefig(os.path.join(fig_dir, f"{args.celltype}_cosine_similarity_heatmap.png"), dpi=300)
            plt.close()
    except Exception as e:
        print(f"Could not compute cosine similarity heatmap: {e}")

if __name__ == "__main__":
    main()

