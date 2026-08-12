#!/usr/bin/env python3
"""
compare_gene_expr.py

Compare gene-expression gene lists (or detected genes) between Cell2location and RCTD
on a per-sample, per-area, per-celltype basis.

Expect directory layout examples:
  RCTD/
    Spatial_1/
      gene_expr_ct/
        dorsal_gm/
          Spatial_1_Astrocytes_gene_expr.csv
          Spatial_1_MI_gene_expr.csv
  cell2location_map/
    Spatial_1/
      gene_expr_ct_mean/
        dorsal_gm/
          Spatial_1_Astrocytes_gene_expr.csv
"""

import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from scipy.stats import spearmanr

def list_samples(base_dir):
    """List sample directories directly under base_dir (e.g. RCTD/Spatial_1)."""
    if not os.path.isdir(base_dir):
        return []
    entries = [d for d in os.listdir(base_dir)
               if os.path.isdir(os.path.join(base_dir, d))]
    return sorted(entries)

def extract_celltype_from_filename(filename, sample_prefix):
    """
    Given filename like "Spatial_1_Astrocytes_gene_expr.csv" and sample_prefix "Spatial_1",
    return "Astrocytes".
    """
    name = os.path.basename(filename)
    if not name.startswith(sample_prefix + "_"):
        return None
    # remove prefix and suffix
    core = name[len(sample_prefix) + 1:]
    # remove suffixes like "_gene_expr.csv" (be permissive)
    core = core.replace("_gene_expr.csv", "")
    core = core.replace(".csv", "")
    return core

def read_genes_from_expr_csv(path):
    """
    Read a gene-expression CSV and return set of genes.
    We assume genes are either in the index or in a column named 'gene'/'Gene'.
    """
    try:
        df = pd.read_csv(path, index_col=0)
        # if index looks like genes, return index
        genes = set(map(str, df.index.astype(str)))
        if len(genes) > 0:
            return genes
    except Exception:
        pass

    # fallback: try a column named gene / Gene
    try:
        df = pd.read_csv(path)
        for col in df.columns:
            if col.lower() == "gene":
                return set(df[col].astype(str))
    except Exception:
        pass

    # if everything fails, return empty set
    return set()

def read_expr_dataframe(path):
    """Return dataframe with genes as index and expression values as a column."""
    try:
        df = pd.read_csv(path, index_col=0)
        # If multiple columns, pick first numeric column
        numeric_cols = df.select_dtypes(include='number').columns
        if len(numeric_cols) > 0:
            return df[[numeric_cols[0]]].rename(columns={numeric_cols[0]: "expression"})
        else:
            return df
    except Exception:
        try:
            df = pd.read_csv(path)
            gene_col = None
            for col in df.columns:
                if col.lower() == "gene":
                    gene_col = col
                    break
            if gene_col is not None:
                df = df.set_index(gene_col)
                numeric_cols = df.select_dtypes(include='number').columns
                if len(numeric_cols) > 0:
                    return df[[numeric_cols[0]]].rename(columns={numeric_cols[0]: "expression"})
        except Exception:
            return pd.DataFrame()
    return pd.DataFrame()


def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)
def plot_rank_scatter(rctd_df, c2l_df, sample, area, ct, out_dir):
    # Keep only shared genes
    common_genes = list(set(rctd_df.index) & set(c2l_df.index))
    if len(common_genes) == 0:
        return
    rctd_ranks = rctd_df.loc[common_genes].sort_values("expression", ascending=False)
    c2l_ranks = c2l_df.loc[common_genes].sort_values("expression", ascending=False)
    ranks_df = pd.DataFrame({
        'gene': common_genes,
        'rctd_rank': [rctd_ranks.index.get_loc(g)+1 for g in common_genes],
        'c2l_rank': [c2l_ranks.index.get_loc(g)+1 for g in common_genes]
    })
    plt.figure(figsize=(4,4))
    plt.scatter(ranks_df['rctd_rank'], ranks_df['c2l_rank'], alpha=0.5)
    max_rank = max(ranks_df['rctd_rank'].max(), ranks_df['c2l_rank'].max())
    plt.plot([0,max_rank],[0,max_rank], 'r--')
    plt.xlabel('RCTD rank')
    plt.ylabel('Cell2location rank')
    plt.title(f'{sample} | {area} | {ct}')
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir,f"{sample}_{area}_{ct}_rank_scatter.png"), dpi=300)
    plt.close()

def plot_topN_overlap(rctd_df, c2l_df, sample, area, ct, out_dir, maxN=200):
    common_genes = list(set(rctd_df.index) & set(c2l_df.index))
    if len(common_genes) == 0:
        return
    rctd_ranks = rctd_df.loc[common_genes].sort_values("expression", ascending=False)
    c2l_ranks = c2l_df.loc[common_genes].sort_values("expression", ascending=False)
    topNs = list(range(10, min(maxN, len(common_genes)), 10))
    overlap_fraction = []
    for N in topNs:
        top_rctd = set(rctd_ranks.index[:N])
        top_c2l = set(c2l_ranks.index[:N])
        overlap_fraction.append(len(top_rctd & top_c2l)/N)
    plt.figure(figsize=(4,3))
    plt.plot(topNs, overlap_fraction, marker='o')
    # Spearman correlation
    ranks_df = pd.DataFrame({
        'rctd_rank':[rctd_ranks.index.get_loc(g)+1 for g in common_genes],
        'c2l_rank':[c2l_ranks.index.get_loc(g)+1 for g in common_genes]
    })
    rho,_ = spearmanr(ranks_df['rctd_rank'], ranks_df['c2l_rank'])
    plt.title(f'{sample} | {area} | {ct}\nSpearman rho={rho:.2f}')
    plt.xlabel('Top N genes')
    plt.ylabel('Fraction overlap')
    plt.ylim(0,1)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir,f"{sample}_{area}_{ct}_topN_overlap.png"), dpi=300)
    plt.close()


def main(args):
    area = args.area
    cell2loc_base = args.cell2loc_dir
    rctd_base = args.rctd_dir
    out_base = args.output_dir

    if not os.path.isdir(rctd_base):
        print(f"[ERROR] RCTD base directory does not exist: {rctd_base}")
        return
    if not os.path.isdir(cell2loc_base):
        print(f"[WARN] cell2location base directory does not exist: {cell2loc_base} (will skip missing samples)")
    ensure_dir(out_base)

    samples = list_samples(rctd_base)
    if len(samples) == 0:
        print(f"[WARN] No sample subdirectories found under RCTD base: {rctd_base}")
        return

    for sample in samples:
        rctd_area_dir = os.path.join(rctd_base, sample, "gene_expr_ct", area)
        c2l_area_dir = os.path.join(cell2loc_base, sample, "gene_expr_ct_mean", area)

        if not os.path.isdir(rctd_area_dir):
            print(f"[WARN] RCTD area directory missing for sample {sample}, area {area}: {rctd_area_dir}")
            continue
        if not os.path.isdir(c2l_area_dir):
            print(f"[WARN] cell2location area directory missing for sample {sample}, area {area}: {c2l_area_dir}")
            # We'll still attempt to compute intersections (which will be empty) but skip detailed plots
            # continue

        # List files in RCTD area dir and extract celltypes
        rctd_files = [f for f in os.listdir(rctd_area_dir) if f.endswith(".csv")]
        rctd_celltypes = {}
        for f in rctd_files:
            ct = extract_celltype_from_filename(f, sample)
            if ct:
                rctd_celltypes[ct] = os.path.join(rctd_area_dir, f)

        if len(rctd_celltypes) == 0:
            print(f"[WARN] No RCTD gene_expr files found for {sample} / {area}")
            continue

        # Determine celltypes present in cell2location for same sample/area (we restrict comparisons to RCTD celltypes)
        c2l_files = {}
        if os.path.isdir(c2l_area_dir):
            for f in os.listdir(c2l_area_dir):
                if f.endswith(".csv"):
                    ct = extract_celltype_from_filename(f, sample)
                    if ct:
                        c2l_files[ct] = os.path.join(c2l_area_dir, f)

        # For this run we only compare the celltypes present in RCTD (CSIDE), and only if corresponding file exists in c2l
        common_celltypes = sorted([ct for ct in rctd_celltypes.keys() if ct in c2l_files])

        sample_out_dir = os.path.join(out_base, sample)
        ensure_dir(sample_out_dir)

        summary_rows = []
        if len(common_celltypes) == 0:
            print(f"[WARN] No common celltypes for sample {sample}, area {area} (RCTD celltypes: {list(rctd_celltypes.keys())}; cell2loc celltypes: {list(c2l_files.keys())})")
            # still write an empty summary csv
            pd_summary = pd.DataFrame(columns=["celltype","n_rctd_genes","n_c2l_genes","n_shared","jaccard"])
            pd_summary.to_csv(os.path.join(sample_out_dir, "overlap_summary.csv"), index=False)
            continue

        for ct in common_celltypes:
            rctd_path = rctd_celltypes[ct]
            c2l_path = c2l_files[ct]

            genes_rctd = read_genes_from_expr_csv(rctd_path)
            genes_c2l  = read_genes_from_expr_csv(c2l_path)

            # Venn plot
            plt.figure(figsize=(4,4))
            venn2([genes_rctd, genes_c2l], set_labels=('RCTD', 'Cell2location'))
            plt.title(f"{sample} | {area} | {ct}")
            out_png = os.path.join(sample_out_dir, f"{sample}_{area}_{ct}_geneexpr_venn.png")
            plt.savefig(out_png, dpi=300, bbox_inches='tight')
            plt.close()

            # Read full expression for rank comparison
            rctd_df = read_expr_dataframe(rctd_path)
            c2l_df = read_expr_dataframe(c2l_path)
            if not rctd_df.empty and not c2l_df.empty:
                plot_rank_scatter(rctd_df, c2l_df, sample, area, ct, sample_out_dir)
                plot_topN_overlap(rctd_df, c2l_df, sample, area, ct, sample_out_dir)

            n_r = len(genes_rctd)
            n_c = len(genes_c2l)
            n_shared = len(genes_rctd & genes_c2l)
            jacc = n_shared / float(len(genes_rctd | genes_c2l)) if (len(genes_rctd | genes_c2l) > 0) else 0.0
            summary_rows.append({
                "celltype": ct,
                "n_rctd_genes": n_r,
                "n_c2l_genes": n_c,
                "n_shared": n_shared,
                "jaccard": jacc,
                "rctd_path": rctd_path,
                "c2l_path": c2l_path
            })

        # Save summary CSV for this sample
        import pandas as pd
        pd.DataFrame(summary_rows).to_csv(os.path.join(sample_out_dir, "overlap_summary.csv"), index=False)
        print(f"[INFO] Wrote {len(summary_rows)} venn plots + summary for sample {sample}, area {area} -> {sample_out_dir}")

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Compare gene_expr per sample/area between Cell2location and RCTD (CSIDE).")
    p.add_argument("--area", required=True, help="Area name (e.g. dorsal_gm)")
    p.add_argument("--cell2loc_dir", required=True, help="Base dir for cell2location samples (e.g. RUN_BASE/cell2location_map)")
    p.add_argument("--rctd_dir", required=True, help="Base dir for RCTD/CSIDE samples (e.g. RUN_BASE/RCTD)")
    p.add_argument("--output_dir", required=True, help="Output directory (plots written to output_dir/<sample>/...)")
    p.add_argument("--ref_label", default=None, help="(optional) ref label")
    args = p.parse_args()
    main(args)

