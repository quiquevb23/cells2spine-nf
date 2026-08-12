'''
    Python file to convert h5ad spatial ref to input required for RCTD
'''

import os
import argparse
import scanpy as sc
import scipy.sparse as sp
import pandas as pd
import multiprocessing
import re
from scipy.io import mmwrite

def parse_args():
    parser = argparse.ArgumentParser(description="Convert scRNA-seq reference for RCTD")
    parser.add_argument("--sc_ref_path", required=True, help="Path to the input h5ad single-cell reference file")
    parser.add_argument("--output_base_dir", required=True, help="Output for files needed for RCTD")
    parser.add_argument('--ref_label', type=str, required=True)
    parser.add_argument("--conditions", help="Comma-separated list of conditions")
    parser.add_argument('--masked_celltypes', type=str, default="")

    return parser.parse_args()

def process_h5ad(adata, label, reference_out_dir):
    print(f"Processing reference: {label}")

    adata.X = adata.layers['counts'].copy()
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    counts = adata.X.T

    mmwrite(os.path.join(reference_out_dir, "counts.mtx"), counts)
    pd.DataFrame(adata.var_names).to_csv(os.path.join(reference_out_dir, "genes.csv"), index=False, header=False)
    pd.DataFrame(adata.obs_names).to_csv(os.path.join(reference_out_dir, "cells.csv"), index=False, header=False)
    
    if ref_label in adata.obs:
        adata.obs[ref_label] = adata.obs[ref_label].str.replace('/', '-', regex=False)
        adata.obs[[ref_label]].to_csv(os.path.join(reference_out_dir, "cell_types.csv"), header=False)

def sanitize_label(label):
    return re.sub(r'[\\/:*?"<>|]', '_', label)

if __name__ == "__main__":
    args = parse_args()
    ref_label = args.ref_label
    out_dir = os.path.join(args.output_base_dir, "data")
    os.makedirs(out_dir, exist_ok=True)

    conditions = [c.strip() for c in args.conditions.split(",") if c.strip()]
    masked_celltypes_raw = [x.strip() for x in args.masked_celltypes.split(",") if x.strip()]
    masked_celltypes = [sanitize_label(ct) for ct in masked_celltypes_raw]

    adata_ref = sc.read_h5ad(args.sc_ref_path)

    for cond in conditions:
        print(f"Found conditions: {conditions}")
        adata_ref_cond = adata_ref[adata_ref.obs['condition'] == cond].copy()
        adata_ref_cond = adata_ref_cond[~adata_ref_cond.obs[ref_label].isin(masked_celltypes)].copy()

        reference_out_dir = os.path.join(out_dir, f"reference_{cond}")
        os.makedirs(reference_out_dir, exist_ok=True)

        process_h5ad(adata_ref_cond, ref_label, reference_out_dir)
    


    #This should be changed to 'condition' or 'field usd for comparisons'
    '''
    subsets = {
        "injured": adata[adata.obs['sample'] == 'Sample_3'].copy(),
        "healthy": adata[adata.obs['sample'] == 'Sample_1'].copy()
    }

    for label, sub_adata in subsets.items():
        reference_out_dir = os.path.join(out_dir, f"reference_{label}")
        os.makedirs(reference_out_dir, exist_ok=True)

        process_h5ad(sub_adata, label, reference_out_dir)
    '''

