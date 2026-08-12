'''
    Python file to convert h5ad spatial h5ad files to input required for RCTD
'''

import os
import argparse
import scanpy as sc
import scipy.sparse as sp
import pandas as pd
import multiprocessing
from scipy.io import mmwrite


def parse_args():
    parser = argparse.ArgumentParser(description="Convert spatial h5ad to RCTD format")
    parser.add_argument("--spatial_input_dir", required=True, help="Directory with sp data after deconvolution with cell2location, input data")
    parser.add_argument("--tmp_dir", required=True, help="Temporary output directory for CSIDE input")
    parser.add_argument(
        '--manual_celltypes', type=str, nargs='+', default=[],
        help="List of cell types to use if --celltype_key is empty. Provide as space-separated list."
    )
    parser.add_argument(
        '--samples', type=str, nargs='+', default=[],
        help="List of sample names."
    )

    return parser.parse_args()

def process(file_path, tmp_dir, sample_name, celltypes):
    adata = sc.read_h5ad(file_path)

    counts = adata.layers['counts'].T if 'counts' in adata.layers else adata.X.T
    counts = sp.csr_matrix(counts)
    df = pd.DataFrame.sparse.from_spmatrix(counts, index=adata.var_names, columns=adata.obs_names)
    df.to_csv(os.path.join(tmp_dir, f"{sample_name}_counts.csv"))

    # Save coordinates for each sample
    if 'array_row' in adata.obs and 'array_col' in adata.obs:
        adata.obs[['array_row', 'array_col']].to_csv(os.path.join(tmp_dir, f"{sample_name}_coordinates.csv"))

    # Extract cell type weights
    if "q05_cell_abundance_w_sf" in adata.obsm:
        weights = pd.DataFrame(
            adata.obsm["q05_cell_abundance_w_sf"],
            index=adata.obs_names,
            columns=adata.obsm["q05_cell_abundance_w_sf"].columns if hasattr(adata.obsm["q05_cell_abundance_w_sf"], "columns") else None

        )
        weights.to_csv(os.path.join(tmp_dir, f"{sample_name}_weights.csv"))
    
    # Extract mapping of 'manual_delineation'
    if "manual_delineation" in adata.obs:
        manual_delineation = pd.DataFrame(
            adata.obs["manual_delineation"],
            index=adata.obs_names,
            columns=adata.obs["manual_delineation"].columns if hasattr(adata.obs["manual_delineation"], "columns") else None
        )
        manual_delineation.to_csv(os.path.join(tmp_dir, f"{sample_name}_manual_delineation.csv"))


    '''
    #Per celltype gene matrices not needed could be needed in other programs (cellcell comm)
    for ct in celltypes:
        counts = adata.obsm[ct].T
        counts = sp.csr_matrix(counts)
        df = pd.DataFrame.sparse.from_spmatrix(counts, index=adata.var_names, columns=adata.obs_names)
        df.to_csv(os.path.join(tmp_dir, f"{sample_name}_{ct}_counts.csv"))
    '''

if __name__ == "__main__":
    args = parse_args()
    os.makedirs(args.tmp_dir, exist_ok=True)

    sample_dir = args.spatial_input_dir
    samples = [i for i in os.listdir(sample_dir) if i.endswith('h5ad')]
    for sample in samples:
        sample_path = os.path.join(sample_dir, sample)
        sample_name = sample.removeprefix("sp").removesuffix(".h5ad")
        print(f"Processing: {sample_path} as sample {sample_name}")
        process(sample_path, args.tmp_dir, sample_name, args.manual_celltypes)

