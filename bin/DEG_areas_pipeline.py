"""
DEA by manually delineated areas using pseudobulk edgeR
"""

import os
import json
import argparse
from pathlib import Path

from scipy import sparse

import pandas as pd
import scanpy as sc
import anndata as ad
import numpy as np
import matplotlib.pyplot as plt

from rpy2 import robjects as ro
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter


# -------------------- ARGUMENTS --------------------
def parse_args():
    parser = argparse.ArgumentParser(description="DEA by manually delineated areas")

    parser.add_argument("--ref_level", required=True)
    parser.add_argument("--delineation_dir", type=str, default="NO_DELINEATION")
    parser.add_argument("--output_base_dir", required=True)
    parser.add_argument("--spatial_input", required=True)
    parser.add_argument("--samples", nargs="+", required=True)
    parser.add_argument("--conditions_map", nargs="+", required=True)
    parser.add_argument("--condition_order", nargs=2, required=True,
        help="Condition contrast order: test first, reference second"
    )

    return parser.parse_args()


args = parse_args()

# -------------------- CONDITIONS MAP --------------------
condition_map = {}
for item in args.conditions_map:
    if "=" in item:
        k, v = item.split("=", 1)
    elif ":" in item:
        k, v = item.split(":", 1)
    else:
        raise ValueError(f"Invalid format: {item}. Expected key=value or key:value")
    condition_map[k] = v


# -------------------- DIRECTORIES --------------------
output_base_dir = Path(args.output_base_dir)
data_dir = output_base_dir / "data"

spatial_input = Path(args.spatial_input)

output_dir = output_base_dir / "DEA_areas_good"
output_dir.mkdir(parents=True, exist_ok=True)

# -------------------- LOAD SAMPLES --------------------
adatas = []

for sample in args.samples:
    print(f"Processing sample: {sample}")

    # Flat h5ad, matching conversor.py / cell2loc_owndata.py's layout
    adata_path = spatial_input / f"{sample}.h5ad"
    if not adata_path.exists():
        raise FileNotFoundError(adata_path)
    adata = sc.read_h5ad(adata_path)

    if args.delineation_dir != "NO_DELINEATION":
        delin_file = data_dir_placeholder = Path(args.delineation_dir) / f"{sample}_manual_delineation.csv"
        if not delin_file.exists():
            raise FileNotFoundError(delin_file)
        df_delineation = pd.read_csv(delin_file, index_col=0)
        adata.obs["manual_delineation"] = df_delineation.reindex(adata.obs_names).iloc[:, 0]
    else:
        adata.obs["manual_delineation"] = "WHOLE_SAMPLE"

    adatas.append(adata)


# -------------------- DEA FUNCTION --------------------
def do_DEA(combined_adata, area_output_dir, area):
    print(f"\n[DEBUG] Entering do_DEA for area {area}")
    print("  Combined shape:", combined_adata.shape)
    print("  Total counts:", combined_adata.X.sum())
    print("  Samples:", combined_adata.obs["sample"].value_counts().to_dict())
    print("  Conditions:", combined_adata.obs["condition"].value_counts().to_dict())

    df = combined_adata.to_df(layer="counts")
    if df.sum().sum() == 0:
        print(f"Area {area} has all-zero counts. Skipping.")
        return

    metadata = combined_adata.obs[["sample", "condition"]]

    grouped = df.groupby(metadata["sample"]).sum()
    grouped_t = grouped

    print("\n[DEBUG] After pseudobulk")
    print(grouped_t.iloc[:5, :5])
    print("Library sizes:", grouped_t.sum(axis=1))

    # ---------------- REMOVE ZERO-LIBRARY SAMPLES ----------------
    library_sizes = grouped_t.sum(axis=1)
    nonzero_samples = library_sizes[library_sizes > 0].index

    print("[DEBUG] Nonzero samples:", list(nonzero_samples))
    print("[DEBUG] All samples:", list(grouped_t.index))
    if len(nonzero_samples) < 2:
        print(f"Skipping area {area}: fewer than 2 samples with non-zero counts")
        return

    grouped_t = grouped_t.loc[nonzero_samples]
    metadata = metadata.loc[metadata["sample"].isin(nonzero_samples)]

    pseudobulk_metadata = pd.DataFrame(grouped_t.index, columns=["sample"])
    pseudobulk_metadata["condition"] = (
        metadata.groupby("sample")["condition"].first().values
    )

    pseudobulk_adata = ad.AnnData(
        X=grouped_t.values,
        obs=pseudobulk_metadata,
        var=pd.DataFrame(index=grouped.columns),
    )

    counts = pd.DataFrame(
        pseudobulk_adata.X,
        index=pseudobulk_adata.obs.index,
        columns=pseudobulk_adata.var_names,
    ).T

    pandas2ri.activate()
    edgeR = importr("edgeR")

    with localconverter(ro.default_converter + pandas2ri.converter):
        r_counts = ro.conversion.py2rpy(counts)
        r_metadata = ro.conversion.py2rpy(pseudobulk_metadata)
        r_genes = ro.conversion.py2rpy(pseudobulk_adata.var_names)


    # ---------------- FORCE CONDITION ORDER ----------------
    test_cond = args.condition_order[0]   # test (numerator)
    ref_cond  = args.condition_order[1]   # reference (denominator)

    r.assign("r_metadata", r_metadata)

    r(f"""
    r_metadata$condition <- factor(
        r_metadata$condition,
        levels = c("{ref_cond}", "{test_cond}")
    )
    """)

    # sanity check (optional but very useful)
    r(f'print(levels(r_metadata$condition))')

    r_metadata = r("r_metadata")

    # Keep r_metadata as R object
    r_condition = r("r_metadata$condition")

    #r_condition = r_metadata.rx2("condition")

    design = r("model.matrix")(ro.Formula("~ condition"), data=r_metadata)

    dge = edgeR.DGEList(counts=r_counts, group=r_condition, genes=r_genes)
    keep = edgeR.filterByExpr(dge, group=r_condition)
    r.assign("dge", dge)
    r.assign("keep", keep)

    dge = r("dge[keep, , keep.lib.sizes=FALSE]")
    dge = edgeR.normLibSizes(dge)
    dge = edgeR.estimateDisp(dge)

    de = edgeR.exactTest(dge)
    results = edgeR.topTags(de, n=r("nrow")(dge)).rx2("table")

    with localconverter(ro.default_converter + pandas2ri.converter):
        results_df = ro.conversion.rpy2py(results)

    # Standardize column names for downstream compatibility
    rename_map = {}

    if "genes" in results_df.columns:
        rename_map["genes"] = "gene"

    if "FDR" in results_df.columns:
        rename_map["FDR"] = "adj.P.Val"

    results_df.rename(columns=rename_map, inplace=True)

    area_output_dir.mkdir(exist_ok=True)
    results_df.to_csv(area_output_dir / f"DEA_{area}_pseudobulk.csv", index=False)

    print(f"Saved DEA results for area {area}")


# -------------------- RUN PER AREA --------------------
all_areas = set()
for adata in adatas:
    all_areas.update(adata.obs["manual_delineation"].unique())
    print(all_areas)

for area in all_areas:
    area_adatas = []
    conditions_present = set()

    for adata in adatas:
        subset = adata[adata.obs["manual_delineation"].astype(str) == str(area)].copy()
        if subset.n_obs == 0:
            continue  # skip empty sample

        # Use raw counts if available
        if "counts" in subset.layers:
            counts_matrix = subset.layers["counts"]
        else:
            counts_matrix = subset.X

        # Convert sparse → dense
        if sparse.issparse(counts_matrix):
            counts_matrix = counts_matrix.toarray()

        # Ensure ndarray (avoid np.matrix)
        counts_matrix = np.asarray(counts_matrix)

        if counts_matrix.ndim != 2:
            print(f"Skipping malformed counts matrix for area {area}, shape={counts_matrix.shape}")
            continue

        # Remove zero-count spots
        spot_sums = counts_matrix.sum(axis=1)
        nonzero_mask = spot_sums > 0
        
        if not np.any(nonzero_mask):
            print(f"[DEBUG] All zero spots for area {area} in sample {subset.obs['sample'].iloc[0]}")
            continue  # skip if no non-zero spots

        subset = subset[nonzero_mask].copy()
        # IMPORTANT: update BOTH X and layer
        subset.layers["counts"] = counts_matrix[nonzero_mask, :]
        subset.X = subset.layers["counts"]


        print(f"[DEBUG] Area {area}, sample {subset.obs['sample'].iloc[0]}")
        print("  shape:", subset.shape)
        print("  total counts:", subset.X.sum())
        print("  nonzero spots:", subset.n_obs)


        area_adatas.append(subset)
        conditions_present.add(subset.obs["condition"].unique()[0])

    # Only run DEA if this area exists in more than 1 condition
    if len(conditions_present) < 2:
        print(f"Skipping area {area}: not represented in both conditions")
        continue

    combined = ad.concat(
        area_adatas,
        join="outer",
        fill_value=0,
    )

    combined.layers["counts"] = combined.X.copy()
    combined.obs_names_make_unique()

    print(f"\n[DEBUG] Combined area {area}")
    print("  Combined shape:", combined.shape)
    print("  Total counts:", combined.X.sum())
    print("  Samples:", combined.obs["sample"].value_counts().to_dict())
    print("  Conditions:", combined.obs["condition"].value_counts().to_dict())

    area_dir = output_dir / str(area)
    area_dir.mkdir(exist_ok=True)

    do_DEA(combined, area_dir, area)

