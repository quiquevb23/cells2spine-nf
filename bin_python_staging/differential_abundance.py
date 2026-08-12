#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat
from sccoda.model.other_models import SimpleModel
import arviz as az
import re

# ---------------------------
# Helper functions
# ---------------------------

def get_present_celltypes(base_dir, sample, area):
    """Return a list of celltypes present for a given sample and area."""
    path = base_dir
    if not os.path.exists(path):
        return []
    celltypes = []
    for f in os.listdir(path):
        if f.endswith(".csv"):
            # Match filenames like: Spatial_1_Astrocytes_gene_expr.csv
            # Remove sample prefix and trailing _gene_expr*.csv
            name = re.sub(rf"^{re.escape(sample)}_", "", f)
            name = re.sub(r"_gene_expr.*\.csv$", "", name)
            celltypes.append(name)
    return celltypes

def renormalize_abundance(df):
    """
    Renormalize 'abundance' per spot so the remaining cell types sum to 1.
    """
    df = df.copy()
    df["abundance"] = df.groupby("spot")["abundance"].transform(lambda x: x / x.sum())
    return df

def load_cell2loc(weights_file, sample):
    df = pd.read_csv(weights_file, index_col=0)
    # No need for normalization of weights
    df = df.reset_index().melt(id_vars="index", var_name="cell_type", value_name="abundance")
    df = df.rename(columns={"index": "spot"})
    df["sample"] = sample
    df["method"] = "cell2loc"
    return df

def load_rctd(weights_file, sample, min_weight=1e-3):
    df = pd.read_csv(weights_file, index_col=0)
    # assume already normalized (RCTD gives already normalized values)
    # Zero out background floor
    df[df < min_weight] = 0.0

    df = df.reset_index().melt(id_vars="index", var_name="cell_type", value_name="abundance")
    df = df.rename(columns={"index": "spot"})
    df["sample"] = sample
    df["method"] = "RCTD"
    return df

def load_area(area_file):
    df = pd.read_csv(area_file)
    df = df.rename(columns={"manual_delineation": "area"})
    df = df.rename(columns={df.columns[0]: "spot"})
    return df

def run_scCODA(df, covariate_col, ref_celltype, output_prefix, test_cond, ref_cond):
    try:
        df[covariate_col] = df[covariate_col].astype(str)

        # Pivot to wide format: rows = spots, columns = cell types
        #comp_df = df.pivot_table(index="spot", columns="cell_type", values="abundance", fill_value=0)
        #comp_df = df.set_index([df.index if "spot" not in df.columns else "spot"])  # index already sample × condition
        #comp_df = comp_df.drop(columns=[covariate_col], errors="ignore")

        comp_df = df.drop(columns=[covariate_col])
        # Make sure counts are numeric
        comp_df = comp_df.apply(pd.to_numeric, errors="coerce").fillna(0)

        # Extract covariates per spot
        #cov_df = df.drop_duplicates(subset="spot")[[ "spot", covariate_col, "total_abundance" ]]
        #cov_df = cov_df.set_index("spot")

        # Extract covariates
        cov_df = df[[covariate_col]].copy()

        # Merge counts and covariate into one DataFrame
        merged_df = comp_df.join(cov_df)

        # 2. Drop the 'sample' column if it exists
        if 'sample' in merged_df.columns:
            merged_df = merged_df.drop(columns=['sample'])

        merged_df[covariate_col] = pd.Categorical(
            merged_df[covariate_col],
            categories=[ref_cond, test_cond],
            ordered=True
        )

        # Ensure the covariate column is categorical string
        #merged_df[covariate_col] = merged_df[covariate_col].astype(str)

        # ✅ Ensure numeric columns first, then covariate last
        import numpy as np
        celltype_cols = merged_df.select_dtypes(include=[np.number]).columns.tolist()
        merged_df = merged_df[celltype_cols + [covariate_col]]

        # Multiply by a scaling factor to represent 'pseudo-cells'
        # This helps the Bayesian model understand the precision of your estimates

        print(f"🧩 Data going into scCODA for {output_prefix}:")
        print(merged_df.dtypes)
        print(merged_df)
        print("Rows (samples):", merged_df.shape[0])

        print("Unique covariate values:", merged_df[covariate_col].unique())
        print("Covariate types:", set(type(x) for x in merged_df[covariate_col]))

        # Create the scCODA dataset (this is the correct syntax for your version)
        data_cc = dat.from_pandas(merged_df, covariate_columns=[covariate_col])
        
        # Fit the model (CHANGE)
        model = mod.CompositionalAnalysis(
            data_cc,
            formula=covariate_col,
            reference_cell_type=ref_celltype
        )
        
        '''
        # 1. Define Cell Types (all columns except 'sample' and the covariate column)
        exclude_cols = ['sample', covariate_col]
        cell_types = [c for c in merged_df.columns if c not in exclude_cols]

        # 2. Extract Data Matrix (Counts)
        # Note: scCODA models usually expect integer counts. 
        # If your floats are whole numbers, we cast to float32 as in your example.
        data_matrix = merged_df[cell_types].values.astype("float32")

        # 3. Build Covariate Matrix
        # formula should be something like "covariate" (no tilde needed for dmatrix sometimes, 
        # but "~ covariate" is the standard R-style Patsy expects)
        formula_str = f"~ {covariate_col}"
        full_cov_matrix = pt.dmatrix(formula_str, merged_df, return_type='dataframe')

        # 4. Handle Reference/Intercept
        # Your example removes the intercept (column 0) to get just the effects
        covariate_names = full_cov_matrix.columns[1:].tolist()
        covariate_matrix_final = full_cov_matrix.iloc[:, 1:].values

        # 5. Define Reference Index
        # If you want the last cell type as the reference:
        #ref_idx = len(cell_types) - 1
        ref_idx = cell_types.index('MOL')
        print(f"Cell types: {cell_types}")
        print(f"Data matrix shape: {data_matrix.shape}")
        print(f"Covariate matrix shape: {covariate_matrix_final.shape}")
        print(f"Reference Index: {ref_idx} (Name: {cell_types[ref_idx]})")    

        # Fit Dirichlet–Multinomial model
        # -----------------------------
        model = SimpleModel(
            covariate_matrix=covariate_matrix_final, 
            data_matrix=data_matrix,
            cell_types=cell_types, 
            covariate_names=covariate_names, 
            formula=formula_str,
            reference_cell_type=ref_idx
        )
        '''
        # Run sampling
        #result = mod.sample_hmc()
        result = model.sample_hmc()
        # Save pickle object
        save_path = f"{output_prefix}_scCODA_result.pkl"
        result.save(save_path)

        # Save summary
        summary = result.summary()
        if summary is not None:
            summary.to_csv(f"{output_prefix}_scCODA_summary.csv")
            return summary
        else:
            print(f"No significant results for {output_prefix} — skipping save.")
            return None

    except Exception as e:
        print(f"scCODA failed for {output_prefix}: {e}")
        return None

def plot_barplots(df, sample, method, outdir):
    """Stacked barplots of top 5 celltypes per area"""
    for area, sub in df.groupby("area"):
        pivot = sub.groupby(["sample", "cell_type"])["abundance"].mean().unstack(fill_value=0)
        # Keep top 5 celltypes by mean abundance
        top5 = pivot.mean().sort_values(ascending=False).head(5).index
        pivot = pivot[top5]

        ax = pivot.plot(kind="bar", stacked=True, figsize=(8,6), colormap="tab20")
        plt.title(f"{sample} - {method} - {area}")
        plt.ylabel("Mean abundance")
        plt.xlabel("Sample")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"{sample}_{method}_{area}_barplot.png"))
        plt.close()

def compare_scCODA_results(cell2loc_summary, rctd_summary, output_prefix):
    """Compare differential abundance results between cell2loc and RCTD."""
    try:
        c2l = cell2loc_summary.copy()
        rctd = rctd_summary.copy()

        # Both summaries include parameter names like "alpha[celltype]"
        c2l = c2l.reset_index().rename(columns={"index": "param"})
        rctd = rctd.reset_index().rename(columns={"index": "param"})

        # Keep only celltype parameters (ignore intercepts etc.)
        c2l = c2l[c2l["param"].str.contains("alpha")]
        rctd = rctd[rctd["param"].str.contains("alpha")]

        # Extract celltype names
        c2l["celltype"] = c2l["param"].str.extract(r"alpha\[(.*)\]")
        rctd["celltype"] = rctd["param"].str.extract(r"alpha\[(.*)\]")

        # Merge by celltype
        merged = pd.merge(
            c2l[["celltype", "mean"]],
            rctd[["celltype", "mean"]],
            on="celltype",
            suffixes=("_cell2loc", "_RCTD")
        )

        # Compute correlation
        corr = merged["mean_cell2loc"].corr(merged["mean_RCTD"])

        # Plot scatter
        plt.figure(figsize=(6,6))
        plt.scatter(merged["mean_cell2loc"], merged["mean_RCTD"], s=60, alpha=0.8)
        plt.axhline(0, color="grey", linestyle="--")
        plt.axvline(0, color="grey", linestyle="--")
        plt.title(f"scCODA comparison\nr = {corr:.2f}")
        plt.xlabel("cell2loc log-fold-change (mean)")
        plt.ylabel("RCTD log-fold-change (mean)")
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_scCODA_comparison_scatter.png")
        plt.close()

        # Save merged summary
        merged["corr"] = corr
        merged.to_csv(f"{output_prefix}_scCODA_comparison_table.csv", index=False)

        print(f"✅ Comparison saved: {output_prefix}_scCODA_comparison_table.csv (r={corr:.2f})")

    except Exception as e:
        print(f"❌ Comparison failed for {output_prefix}: {e}")


# ---------------------------
# Main
# ---------------------------

def main():
    parser = argparse.ArgumentParser(description="Differential abundance with scCODA on cell2loc and RCTD outputs.")
    parser.add_argument("--input_cell2loc", required=True)
    parser.add_argument("--input_rctd", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--data_dir", required=True)
    parser.add_argument("--ref_label", required=True)
    parser.add_argument("--masked_celltypes", nargs="*", default=[])
    parser.add_argument("--samples", nargs="+", required=True)
    parser.add_argument("--conditions_map", nargs="+", required=True,
                        help="Mapping sample:condition, e.g. Spatial_1:Control Spatial_2:Disease")
    parser.add_argument("--condition_order", nargs=2, required=True,
        help="Condition contrast order: test first, reference second"
    )

    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # Parse condition map
    cond_map = dict(s.split(":") for s in args.conditions_map)

    # Collect per-sample data here, then concatenate to run scCODA across conditions
    all_cell2loc = []
    all_rctd = []
    all_areas = []

    test_cond = args.condition_order[0]   # test (numerator)
    ref_cond  = args.condition_order[1]   # reference (denominator)

    for sample in args.samples:
        cond = cond_map.get(sample, "Unknown")

        # Input files
        cell2loc_file = os.path.join(args.input_cell2loc, sample, f"{sample}_cell2loc_weights.csv")
        rctd_file     = os.path.join(args.input_rctd, sample, f"{sample}_RCTD_weights.csv")
        area_file     = os.path.join(args.data_dir, sample, f"{sample}_manual_delineation.csv")

        if not (os.path.exists(cell2loc_file) and os.path.exists(rctd_file) and os.path.exists(area_file)):
            print(f"Skipping {sample}, missing files")
            continue

        # Load data
        df_cell2loc = load_cell2loc(cell2loc_file, sample)
        df_rctd = load_rctd(rctd_file, sample)
        df_area = load_area(area_file)

        # Merge area info into both dataframes
        df_cell2loc = df_cell2loc.merge(df_area, on="spot", how="left")
        df_rctd     = df_rctd.merge(df_area, on="spot", how="left")

        # Add condition column
        df_cell2loc["condition"] = cond
        df_rctd["condition"] = cond

        # Append to global lists
        all_cell2loc.append(df_cell2loc)
        all_rctd.append(df_rctd)
        all_areas.append(df_area)

    # If nothing loaded, exit
    if len(all_cell2loc) == 0 and len(all_rctd) == 0:
        print("No samples loaded — exiting.")
        return

    # Concatenate across samples
    df_cell2loc_all = pd.concat(all_cell2loc, ignore_index=True) if len(all_cell2loc) else pd.DataFrame()
    df_rctd_all     = pd.concat(all_rctd, ignore_index=True)     if len(all_rctd) else pd.DataFrame()
    df_all_methods  = pd.concat([df_cell2loc_all, df_rctd_all], ignore_index=True)

    # Combined area set (union of all area files)
    if len(all_areas):
        df_area_all = pd.concat(all_areas, ignore_index=True).drop_duplicates().reset_index(drop=True)
    else:
        df_area_all = pd.DataFrame(columns=["spot", "area"])

    # Save merged overview
    df_all_methods.to_csv(os.path.join(args.out_dir, "all_samples_merged_abundances.csv"), index=False)

    # Run scCODA per method × area (comparing conditions across samples)
    summaries = {}
    for method, df_method in [("cell2loc", df_cell2loc_all), ("RCTD", df_rctd_all)]:
        if df_method.empty:
            continue

        # Ensure area column exists; if not, attempt to merge using df_area_all
        if "area" not in df_method.columns and not df_area_all.empty:
            df_method = df_method.merge(df_area_all, on="spot", how="left")

        # iterate areas found in the combined area table OR in the method df
        areas = sorted(set(df_area_all["area"].dropna().unique().tolist() + df_method["area"].dropna().unique().tolist()))
        areas = ['Ventral_gm']
        for area in areas:
            sub = df_method[df_method["area"] == area].copy()

            # Get list of valid celltypes for this area
            base_dir = args.input_cell2loc if method == "cell2loc" else args.input_rctd
            celltype_dir = "gene_expr_ct_mean" if method == "cell2loc" else "gene_expr_ct"
            present_cts = get_present_celltypes(os.path.join(base_dir, sample, celltype_dir, area), sample, area)

            if present_cts:
                sub = sub[sub["cell_type"].isin(present_cts)]

            if args.masked_celltypes:
                sub = sub[~sub["cell_type"].isin(args.masked_celltypes)]

            #sub = renormalize_abundance(sub)
            
            '''
            sub["total_abundance"] = sub.groupby("spot")["abundance"].transform("sum")
            if sub.empty:
                # nothing to run
                continue
            '''
            celltype_cols = sub["cell_type"].unique()  # get all present cell types

            # Pivot back to wide format: rows = spot, columns = cell type
            sub_wide = sub.pivot_table(
                index=["sample", "condition"], 
                columns="cell_type", 
                values="abundance", 
                aggfunc="sum",   # sum over spots in the same sample × condition
                fill_value=0
            ).reset_index()

            if sub_wide.empty:
                continue

            # Need at least two conditions to perform differential testing
            unique_conditions = sub_wide["condition"].dropna().unique()
            if len(unique_conditions) < 2:
                print(f"Skipping scCODA for {method} - {area}: found conditions = {unique_conditions}")
                continue

            out_area_dir = os.path.join(args.out_dir, area)
            os.makedirs(out_area_dir, exist_ok=True)
            out_prefix = os.path.join(out_area_dir, f"{area}_{method}")
            
            summary = run_scCODA(
                df=sub_wide.rename(columns={"condition": "covariate"}),
                covariate_col="covariate",
                ref_celltype="MOL",
                output_prefix=out_prefix,
                test_cond=test_cond,
                ref_cond=ref_cond
                )
            if summary is not None:
                summaries[(method, area)] = summary

        # Plot barplots for the method across samples/areas
        plot_barplots(df_method, "all_samples", method, args.out_dir)

    # Compare scCODA results method-by-method for each area
    all_areas_list = sorted({a for a in df_area_all["area"].dropna().unique()} | {a for a in df_all_methods["area"].dropna().unique()})
    for area in all_areas_list:
        key_c2l = ("cell2loc", area)
        key_rctd = ("RCTD", area)
        if key_c2l in summaries and key_rctd in summaries:
            out_prefix = os.path.join(args.out_dir, area, f"{area}_comparison")
            compare_scCODA_results(summaries[key_c2l], summaries[key_rctd], out_prefix)

    print("✅ Differential abundance (injured vs healthy) + comparison finished.")

if __name__ == "__main__":
    main()

