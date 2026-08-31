#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import numpy as np
from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat

# ---------------------------
# IO helpers
# ---------------------------

def load_cell2loc(weights_file, sample):
    df = pd.read_csv(weights_file, index_col=0)
    df = df.reset_index().melt(
        id_vars="index",
        var_name="cell_type",
        value_name="abundance"
    )
    df = df.rename(columns={"index": "spot"})
    df["sample"] = sample
    return df

def load_rctd(weights_file, sample, min_weight=1e-3):
    df = pd.read_csv(weights_file, index_col=0)
    df[df < min_weight] = 0.0
    df = df.reset_index().melt(
        id_vars="index",
        var_name="cell_type",
        value_name="abundance"
    )
    df = df.rename(columns={"index": "spot"})
    df["sample"] = sample
    return df

def load_area(area_file):
    df = pd.read_csv(area_file)
    df = df.rename(columns={df.columns[0]: "spot",
                            "manual_delineation": "area"})
    return df

# ---------------------------
# scCODA runner
# ---------------------------

def run_scCODA_spotlevel(
    df,
    covariates,
    ref_celltype,
    output_prefix,
    test_cond,
    ref_cond
):
    try:
        # Ensure categorical covariates
        for c in covariates:
            df[c] = df[c].astype(str)

        # Set condition order
        df["condition"] = pd.Categorical(
            df["condition"],
            categories=[ref_cond, test_cond],
            ordered=True
        )

        # Identify count columns
        count_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        merged_df = df[count_cols + covariates]

        print(f"\n🧩 scCODA input ({output_prefix})")
        print(merged_df.dtypes)
        print(merged_df)
        print("Rows:", merged_df.shape[0])

        data_cc = dat.from_pandas(
            merged_df,
            covariate_columns=covariates
        )

        model = mod.CompositionalAnalysis(
            data_cc,
            formula=" + ".join(covariates),
            reference_cell_type=ref_celltype
        )

        result = model.sample_hmc(
            num_results=20000,
            num_burnin=5000,
            verbose=True
        )
        # Save pickle object
        save_path = f"{output_prefix}_scCODA_result_spotlevel.pkl"
        result.save(save_path)


        summary = result.summary()
        if summary is None:
            print(f"⚠️ No credible effects for {output_prefix}")
            return None

        summary.to_csv(f"{output_prefix}_scCODA_summary.csv")
        return summary

    except Exception as e:
        print(f"❌ scCODA failed for {output_prefix}: {e}")
        return None

# ---------------------------
# Main
# ---------------------------

def main():
    parser = argparse.ArgumentParser()
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

    cond_map = dict(s.split(":") for s in args.conditions_map)
    test_cond, ref_cond = args.condition_order

    all_cell2loc = []
    all_rctd = []

    for sample in args.samples:
        cond = cond_map[sample]

        df_c2l = load_cell2loc(
            os.path.join(args.input_cell2loc, sample, f"{sample}_cell2loc_weights.csv"),
            sample
        )
        df_rctd = load_rctd(
            os.path.join(args.input_rctd, sample, f"{sample}_RCTD_weights.csv"),
            sample
        )
        df_area = load_area(
            os.path.join(args.data_dir, sample, f"{sample}_manual_delineation.csv")
        )

        for df in (df_c2l, df_rctd):
            df["condition"] = cond
            df.merge(df_area, on="spot", how="left")

        all_cell2loc.append(df_c2l.merge(df_area, on="spot"))
        all_rctd.append(df_rctd.merge(df_area, on="spot"))

    for method, df_all in {
        "cell2loc": pd.concat(all_cell2loc),
        "RCTD": pd.concat(all_rctd)
    }.items():

        for area in sorted(df_all["area"].dropna().unique()):
            sub = df_all[df_all["area"] == area]

            # Pivot to spot-level wide
            sub_wide = sub.pivot_table(
                index=["spot", "sample", "condition"],
                columns="cell_type",
                values="abundance",
                fill_value=0
            ).reset_index()

            if sub_wide.shape[0] < 20:
                continue

            covariates = ["condition", "sample"]
            #covariates = "condition"
            out_prefix = os.path.join(
                args.out_dir,
                f"{area}_{method}"
            )

            run_scCODA_spotlevel(
                df=sub_wide,
                covariates=covariates,
                ref_celltype="MOL",
                output_prefix=out_prefix,
                test_cond=test_cond,
                ref_cond=ref_cond
            )

    print("✅ Spot-level scCODA finished.")

if __name__ == "__main__":
    main()

