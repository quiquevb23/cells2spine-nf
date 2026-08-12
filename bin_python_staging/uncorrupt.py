import os
import h5py
import shutil
import argparse
import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix

def parse_args():
    parser = argparse.ArgumentParser(description="Clean corrupted .h5ad layers and export cleaned files")
    parser.add_argument("--base_dir", type=str, required=True)
    parser.add_argument("--sp_deconv_data_dir", type=str, required=True)
    parser.add_argument("--ct_gene_expr_output", type=str, required=True)
    parser.add_argument("--manual_celltypes", nargs='*', default=[])
    parser.add_argument("--samples", nargs='*', default=[])
    return parser.parse_args()

def clean_and_reconstruct_layers(sample, input_dir, output_dir):
    original_file = os.path.join(input_dir, f"sp{sample}.h5ad")
    cleaned_file = os.path.join(output_dir, f"sp{sample}_cleaned.h5ad")

    # Make sure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Copy file
    shutil.copyfile(original_file, cleaned_file)

    # Step 1: Clean bad layers
    with h5py.File(cleaned_file, "a") as f:
        if 'obsm' in f:
            bad_layers = []
            for key in list(f['obsm'].keys()):
                if not isinstance(f['obsm'][key], h5py.Dataset):
                    bad_layers.append(key)
            for layer in bad_layers:
                del f['obsm'][layer]
            print(f"{sample}: Removed layers: {bad_layers}")

    # Step 2: Read cleaned AnnData
    adata = sc.read_h5ad(cleaned_file)

    # Step 3: Reconstruct sparse matrices from original file
    reconstructed_layers = {}
    with h5py.File(original_file, "r") as f:
        if 'X' in f and isinstance(f['X'], h5py.Group):
            indptr = f['X']['indptr'][()]
            indices = f['X']['indices'][()]
            shape = (len(indptr) - 1, indices.max() + 1)
        else:
            raise ValueError(f"Cannot determine shape from {original_file}")

        if 'obsm' in f:
            for layer_name in f['obsm']:
                group = f['obsm'][layer_name]
                if isinstance(group, h5py.Dataset):
                    continue
                try:
                    data = group['data'][()]
                    indices = group['indices'][()]
                    indptr = group['indptr'][()]
                    matrix = csr_matrix((data, indices, indptr), shape=shape)
                    reconstructed_layers[layer_name] = matrix
                except Exception as e:
                    print(f"{sample}: Failed to reconstruct {layer_name}: {e}")

    # Step 4: Insert into .obsm
    for layer_name, matrix in reconstructed_layers.items():
        adata.obsm[layer_name] = matrix.toarray()
        print(f"{sample}: Added layer '{layer_name}' to .obsm")

    # Save cleaned file again
    adata.write(cleaned_file)
    print(f"{sample}: Cleaned file saved to {cleaned_file}")

def main():
    args = parse_args()

    input_dir = args.sp_deconv_data_dir
    output_dir = os.path.join(input_dir, "cleaned")
    samples = args.samples
    manual_celltypes = [ct.strip() for ct in args.manual_celltypes if ct.strip()]

    if not samples:
        # Auto-detect samples
        files = os.listdir(input_dir)
        samples = [f.replace("sp", "").replace(".h5ad", "") for f in files if f.startswith("sp") and f.endswith(".h5ad")]

    print(f"Cleaning samples: {samples}")

    for sample in samples:
        clean_and_reconstruct_layers(sample, input_dir, output_dir)

if __name__ == "__main__":
    main()
