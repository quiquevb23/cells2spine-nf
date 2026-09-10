#!/usr/bin/env python3
#"run_l3_20260521_162030_SCI_10d,healthy_rostral"

import os
import scanpy as sc
import pandas as pd
from cell2location.models import RegressionModel

# Directory where your files are located
base_dir = "/storage/gge/Quique/Cells2SpineData/full_exp/spatial/matrices/Deconvolution_owndata/run_l3_20260521_162030_SCI_10d,healthy_rostral/cell2location_map/reference_signatures_SCI_10d/"

# File paths
adata_path = os.path.join(base_dir, "sc.h5ad")
#model_path = os.path.join(base_dir, "model.pt")

# Load AnnData
adata_ref = sc.read_h5ad(adata_path)

# Load model
mod = RegressionModel.load(base_dir, adata_ref)

# Extract inferred average expression per cluster (cell type signatures)
inf_aver = adata_ref.varm['means_per_cluster_mu_fg'].copy()

# Rename columns to cell type names
inf_aver.columns = adata_ref.uns['mod']['factor_names']

# Save DataFrame
output_path = os.path.join(base_dir, "celltype_signatures.csv")
inf_aver.to_csv(output_path)

print(f"Saved cell type signatures to: {output_path}")

