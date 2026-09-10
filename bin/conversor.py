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

def parse_args():
    parser = argparse.ArgumentParser(description="Convert spatial h5ad to RCTD format")
    parser.add_argument("--spatial_input", required=True, help="Directory with spatial h5ad files")
    parser.add_argument("--output_base_dir", required=True, help="Output for files needed for RCTD")
    parser.add_argument("--samples", nargs="+", help="List of spatial sample names")

    return parser.parse_args()

def clean_ct_name(ct):
    return re.sub(r'^meanscell_abundance_w_sf_', '', ct)

def process(file_path, sample_out_dir, sample_name, sample_path):
    adata = sc.read_h5ad(file_path)

    # Save raw spatial counts
    #counts = adata.layers['counts'].T if 'counts' in adata.layers else adata.X.T
    #counts = sp.csr_matrix(counts)
    #df = pd.DataFrame.sparse.from_spmatrix(counts, index=adata.var_names, columns=adata.obs_names)
    #df.to_csv(os.path.join(sample_out_dir, f"{sample_name}_counts.csv"))

    # Save coordinates for each sample
    #if 'array_row' in adata.obs and 'array_col' in adata.obs:
    #    adata.obs[['array_row', 'array_col']].to_csv(os.path.join(sample_out_dir, f"{sample_name}_coords.csv"))

    if 'imagerow' in adata.obs and 'imagecol' in adata.obs:
        adata.obs[['imagerow', 'imagecol']].to_csv(os.path.join(sample_out_dir, f"{sample_name}_coords.csv"))


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

    spatial_input = args.spatial_input

    samples = args.samples

    for sample_name in samples:
        try:
            print(f"Processing {sample_name}...")
            sample_path = os.path.join(spatial_input, sample_name, "outs", "matrices")
            h5ad_file = f"{sample_name}_manual_delineation.h5ad"
            #sample_path = os.path.join(args.output_base_dir, "cell2location_map", sample_name)
            #h5ad_file = f"sp{sample_name}.h5ad"

            full_path = os.path.join(sample_path, h5ad_file)

            if not os.path.exists(full_path):
                print(f"❌ File not found: {full_path}")
                continue

            sample_out_dir = os.path.join(out_dir, sample_name)
            os.makedirs(sample_out_dir, exist_ok=True)
            process(full_path, sample_out_dir, sample_name, sample_path)

            print(f"✅ Finished {sample_name}")

        except Exception as e:
            print(f"🔥 ERROR in {sample_name}: {e}")
            raise  # rethrow so you SEE the failure
