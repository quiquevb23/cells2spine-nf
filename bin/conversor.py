#!/usr/bin/env python3
'''
    Python file to convert h5ad spatial h5ad files to input required for RCTD
'''

import os
import argparse
import scanpy as sc
import scipy.sparse as sp
import pandas as pd
import re
import glob

def parse_args():
    parser = argparse.ArgumentParser(description="Convert spatial h5ad to RCTD format")
    parser.add_argument("--spatial_input", required=True, help="Directory with spatial h5ad files")
    parser.add_argument("--output_base_dir", required=True, help="Output for files needed for RCTD")
    parser.add_argument("--samples", nargs="+", help="List of spatial sample names")
    parser.add_argument("--spatial_h5ad_glob", default="*_manual_delineation.h5ad",
                        help="Glob pattern for spatial h5ad files within each sample directory")
    return parser.parse_args()

def clean_ct_name(ct):
    return re.sub(r'^meanscell_abundance_w_sf_', '', ct)

def write_coords(adata, coords_dir, sample_name):
    """<coords_dir>/<sample>_coords.csv (barcode, imagerow, imagecol) for CSIDE.R / RCTD."""
    if 'imagerow' in adata.obs and 'imagecol' in adata.obs:
        coords = adata.obs[['imagerow', 'imagecol']].copy()
    elif 'spatial' in adata.obsm:
        # obsm['spatial'] is (x, y) = (imagecol, imagerow) in full-res pixels
        xy = adata.obsm['spatial']
        coords = pd.DataFrame({'imagerow': xy[:, 1], 'imagecol': xy[:, 0]}, index=adata.obs_names)
    elif 'array_row' in adata.obs and 'array_col' in adata.obs:
        coords = adata.obs[['array_row', 'array_col']].rename(
            columns={'array_row': 'imagerow', 'array_col': 'imagecol'})
    else:
        raise KeyError(f"No spatial coordinates for {sample_name}: need obs imagerow/imagecol, "
                       f"obsm['spatial'] or obs array_row/array_col")
    coords.to_csv(os.path.join(coords_dir, f"{sample_name}_coords.csv"))
    print(f"[OK] wrote {sample_name}_coords.csv ({len(coords)} spots)")


def write_counts(adata, counts_dir, sample_name):
    """<counts_dir>/<sample>_counts.csv, genes x spots raw counts (read by run_masigpro.R)."""
    counts = adata.layers['counts'] if 'counts' in adata.layers else adata.X
    counts = sp.csr_matrix(counts).T.tocsr()
    df = pd.DataFrame.sparse.from_spmatrix(counts, index=adata.var_names, columns=adata.obs_names)
    df.to_csv(os.path.join(counts_dir, f"{sample_name}_counts.csv"))
    print(f"[OK] wrote {sample_name}_counts.csv ({counts.shape[0]} genes x {counts.shape[1]} spots)")


def process(file_path, sample_out_dir, sample_name, sample_path, coords_dir, counts_dir):
    adata = sc.read_h5ad(file_path)

    write_coords(adata, coords_dir, sample_name)
    write_counts(adata, counts_dir, sample_name)


    # Extract mapping of 'manual_delineation'
    if "manual_delineation" in adata.obs:
        manual_delineation = adata.obs[["manual_delineation"]]
        manual_delineation.to_csv(os.path.join(sample_out_dir, f"{sample_name}_manual_delineation.csv"))

    # ---- NEW: compute and save HVGs ----
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=5000,
        flavor="seurat_v3",
        inplace=True
    )

    # Subset only HVGs
    hvg_mask = adata.var["highly_variable"]
    hvg_df = pd.DataFrame(
        {
            "highly_variable_rank": adata.var.loc[hvg_mask, "highly_variable_rank"]
        },
        index=adata.var.index[hvg_mask]
    )

    # Save HVGs
    hvg_df.to_csv(os.path.join(sample_out_dir, f"{sample_name}_HVGs.csv"))


if __name__ == "__main__":
    args = parse_args()

    out_dir = os.path.join(args.output_base_dir, "data")
    os.makedirs(out_dir, exist_ok=True)
    # EXTRACT_SPATIAL_INPUTS declares these two as its outputs
    coords_dir = os.path.join(args.output_base_dir, "coords")
    counts_dir = os.path.join(args.output_base_dir, "counts")
    os.makedirs(coords_dir, exist_ok=True)
    os.makedirs(counts_dir, exist_ok=True)

    spatial_input = args.spatial_input

    samples = args.samples

    for sample_name in samples:
        try:
            print(f"Processing {sample_name}...")
            sample_path = os.path.join(spatial_input, sample_name, "outs", "matrices")
            matches = glob.glob(os.path.join(sample_path, args.spatial_h5ad_glob))

            if not matches:
                print(f"❌ No h5ad matching '{args.spatial_h5ad_glob}' in {sample_path}")
                continue

            full_path = matches[0]

            sample_out_dir = os.path.join(out_dir, sample_name)
            os.makedirs(sample_out_dir, exist_ok=True)
            process(full_path, sample_out_dir, sample_name, sample_path, coords_dir, counts_dir)

            print(f"✅ Finished {sample_name}")

        except Exception as e:
            print(f"🔥 ERROR in {sample_name}: {e}")
            raise  # rethrow so you SEE the failure

    # Fail loudly instead of handing CSIDE an empty coords dir
    missing = [s for s in samples
               if not os.path.exists(os.path.join(coords_dir, f"{s}_coords.csv"))]
    if missing:
        raise SystemExit(f"No coords written for samples: {missing} (h5ad not found?)")
