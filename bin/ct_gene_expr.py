#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import numpy as np
from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
import matplotlib.pyplot as plt

# ---------------- Functions ---------------- #

def parse_args():
    parser = argparse.ArgumentParser(
        description="Pseudobulk DEA per cell type and delineation area (from CSVs)."
    )
    parser.add_argument(
        "--base_dir", type=str, required=True,
        help="Base run name folder (contains data/ and cell2location_map/)."
    )
    parser.add_argument(
        "--out_dir", type=str, default="CT_Gene_expr",
        help="Output directory for results."
    )
    parser.add_argument(
        "--cside_dir", type=str, default="CSIDE",
        help="Output directory for CSIDE dir."
    )
    parser.add_argument(
        "--samples", type=str, nargs="+", required=True,
        help="List of sample names to include."
    )
    parser.add_argument(
        "--celltypes", type=str, nargs="+", required=True,
        help="List of cell types to process."
    )
    parser.add_argument(
        "--conditions_map", type=str, nargs="+", required=True,
        help="Mapping of samples to conditions, in the form 'sample:condition'."
    )
    parser.add_argument(
        '--pseudocount', type=float, default=1e-3,
        help='Pseudocount added before log2 transform'
    )
    parser.add_argument(
        '--condition_order', nargs=2, required=True,
        help="Condition contrast order: test first, reference second"
    )

    parser.add_argument(
        '--delineation_dir', type=str, default="NO_DELINEATION",
        help="Directory containing <sample>_manual_delineation.csv files. "
             "Use NO_DELINEATION to process all spots as a single WHOLE_SAMPLE region."
    )

    return parser.parse_args()

def parse_conditions_map(conditions_list):
    """
    Convert ['Spatial_1:healthy', 'Spatial_2:disease'] -> {'Spatial_1': 'healthy', 'Spatial_2': 'disease'}
    """
    mapping = {}
    for item in conditions_list:
        if ":" not in item:
            raise ValueError(f"Condition map entry '{item}' is not in 'sample:condition' format")
        sample, condition = item.split(":", 1)
        mapping[sample] = condition
    return mapping

def filter_celltype_area(df_expr, weight_threshold=0.1, min_spot_fraction=0.05,
                         max_weight_threshold=0.2, median_weight_threshold=0.05, verbose=True):
    """
    df_expr: barcodes x genes dataframe (after merging with delineation)
    Returns True if the cell type passes the filter, False otherwise.
    """
    n_spots_total = df_expr.shape[0]
    if n_spots_total == 0:
        if verbose:
            print("[filter] No spots for this cell type.")
        return False

    # Compute summary metrics
    # Using max and median across genes per spot
    max_per_spot = df_expr.max(axis=1)
    median_per_spot = df_expr.median(axis=1)

    # Count number of barcodes above weight_threshold
    n_spots = (max_per_spot >= weight_threshold).sum()

    # Minimum spots threshold
    min_spots_required = int(np.ceil(min_spot_fraction * n_spots_total))

    # Median / max thresholds (heuristics)
    median_val = median_per_spot.median()
    max_val = max_per_spot.max()

    keep = (n_spots >= min_spots_required) and ((max_val >= max_weight_threshold) or (median_val >= median_weight_threshold))

    if verbose:
        print(f"[filter] total spots={n_spots_total}, above threshold={n_spots}, min required={min_spots_required}, max_val={max_val:.3f}, median_val={median_val:.3f}, keep={keep}")
    
    return keep

def run_limma_pseudobulk(counts_df, metadata_df, output_path, pseudocount=1e-3, plot_prefix=None, condition_order=None):
    """
    counts_df: pandas.DataFrame (genes x samples)
    metadata_df: pandas.DataFrame with columns ['sample','condition'] (rows in any order)
    output_path: CSV path to save topTable
    plot_prefix: file path prefix (without extension) where we'll save MA and mean-var plots
    """
    pandas2ri.activate()
    limma = importr('limma')
    base = importr('base')

    # Align counts columns to metadata sample order
    sample_order = list(metadata_df['sample'])
    if set(sample_order) != set(counts_df.columns.tolist()):
        raise ValueError(f"Mismatch between counts samples and metadata samples. counts: {counts_df.columns.tolist()[:10]}..., meta: {sample_order[:10]}...")
    counts_df = counts_df[sample_order]

    #DEBUGS
    #print("[limma] counts_df shape:", counts_df.shape)
    #print("[debug] metadata_df:")
    #print(metadata_df)

    # 1) Replace NA with 0 (absence) and ensure numeric
    counts_df = counts_df.fillna(0.0).astype(float)


    # 2) Remove genes with all zeros across all samples
    nonzero_mask = (counts_df.sum(axis=1) > 0)
    counts_df = counts_df.loc[nonzero_mask, :]
    if counts_df.shape[0] == 0:
        print("[limma] No non-zero genes after filtering. Skipping.")
        return

    # 3) log2 transform
    log_expr = np.log2(counts_df + pseudocount)
    #log_expr = counts_df

    print(log_expr.index[:10])
    print(log_expr.shape)

    # Convert to R objects
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_expr = ro.conversion.py2rpy(log_expr)
        r_meta = ro.conversion.py2rpy(metadata_df)


    # Set R row/col names
    #ro.r.assign("r_expr", r_expr)
    #ro.r.assign("gene_names", list(log_expr.index))
    #ro.r.assign("sample_names", list(log_expr.columns))
    #ro.r("rownames(r_expr) <- gene_names")
    #ro.r("colnames(r_expr) <- sample_names")

    ###NEW BLOCK
    # Assign clean gene & sample names
    gene_names = list(log_expr.index.astype(str))
    sample_names = list(log_expr.columns.astype(str))

    ro.globalenv['r_expr'] = r_expr
    ro.globalenv['gene_names'] = ro.StrVector(gene_names)
    ro.globalenv['sample_names'] = ro.StrVector(sample_names)

    # Properly set row/colnames BEFORE normalization
    ro.r("rownames(r_expr) <- gene_names")
    ro.r("colnames(r_expr) <- sample_names")

    # Now safe to normalize
    ro.r("""
    r_expr_norm <- normalizeBetweenArrays(r_expr, method="quantile")
    rownames(r_expr_norm) <- gene_names       # <- restore real gene names
    colnames(r_expr_norm) <- sample_names
    r_expr_to_use <- r_expr_norm
    """)

    # Correct debugging print
    print(ro.r("head(rownames(r_expr_to_use))"))

    ########

    # 🔧 Ensure condition is factor INSIDE r_meta
    #r_meta = ro.r('transform')(r_meta, condition=ro.r('as.factor')(r_meta.rx2('condition')))
    
    # In R: perform between-array normalization (quantile) ##CHANGED AGAIN:
    #ro.r('r_expr_norm <- normalizeBetweenArrays(r_expr, method="quantile")')
    #ro.r('r_expr_to_use <- r_expr_norm')

    # Normalize but preserve rownames
    #ro.r('''
    #r_expr_norm <- normalizeBetweenArrays(r_expr, method="quantile")
    #rownames(r_expr_norm) <- rownames(r_expr)  # 🔹 preserve gene names
    #r_expr_to_use <- r_expr_norm
    #''')


    # Assign so we can inspect in R
    ro.r.assign("r_meta", r_meta)
    # Ensure condition is a factor IN R
    #ro.r("r_meta$condition <- as.factor(r_meta$condition)")

    #ro.r('r_meta$condition <- relevel(r_meta$condition, ref="healthy")')

    if condition_order is None or len(condition_order) != 2:
        raise ValueError("condition_order must be length 2: [test, reference]")

    test_cond = condition_order[0]
    ref_cond  = condition_order[1]

    ro.globalenv["test_cond"] = test_cond
    ro.globalenv["ref_cond"] = ref_cond

    ro.r("""
    r_meta$condition <- factor(
        r_meta$condition,
        levels = c(ref_cond, test_cond)
    )
    """)

    # sanity check
    ro.r("print(levels(r_meta$condition))")


    print("[debug] r_meta as seen in R:")
    ro.r("print(r_meta)")
    ro.r("print(str(r_meta))")
    
    # put factor back into the metadata object in R for design creation
    #ro.r('condition <- NULL')()
    # build a data.frame in R with sample and condition
    # We'll use the r_meta we already passed to model.matrix; model.matrix(~ condition, data=r_meta)


    # build design in R
    #design = ro.r('model.matrix')(ro.Formula('~ condition'), data=r_meta)
    #design = ro.r("model.matrix(~ condition, data=r_meta)")
    # Build design safely
    ro.r("design <- model.matrix(~ condition, data=r_meta)")

    with localconverter(ro.default_converter + pandas2ri.converter):
        #design_colnames = list(ro.conversion.rpy2py(ro.r('colnames')(design)))
        design_colnames = list(ro.conversion.rpy2py(ro.r("colnames(design)")))
    print('[debug] design colnames:', design_colnames)

    if len(design_colnames) < 2:
        raise ValueError('[limma] Design matrix has <2 columns; cannot determine condition coefficient')

    coef_index = len(design_colnames)  # limma is 1-based in R
    coef_name = design_colnames[-1]
    print(f"[limma] Using coef index {coef_index} (R colname: {coef_name}) for testing")

    # Fit limma
    #fit = limma.lmFit(r_expr, ro.r("design"))
    #fit2 = limma.eBayes(fit, trend=True)
    fit = limma.lmFit(ro.r("r_expr_to_use"), ro.r("design"))
    fit2 = limma.eBayes(fit, trend=True)

    # topTable: request all genes, sort by adj.P.Val
    #n_genes = int(ro.r('nrow')(r_expr)[0])
    n_genes = int(ro.r('nrow(r_expr_to_use)')[0])
    #n_genes = int(ro.r('nrow')(r_expr_to_use)[0])
    
    #top = limma.topTable(fit2, number=n_genes, sort_by='P', adjust_method='BH', coef=coef_index)

    
    # Use a direct R call so genelist works reliably
    ro.globalenv["fit2"] = fit2
    ro.globalenv['n_genes'] = n_genes
    ro.globalenv['coef_index'] = coef_index
    ro.r("""
    results <- topTable(
        fit2,
        number = n_genes,
        sort.by = "P",
        adjust.method = "BH",
        coef = coef_index,
        genelist = rownames(r_expr_to_use)
    )
    """)

    with localconverter(ro.default_converter + pandas2ri.converter):
        res_df = ro.conversion.rpy2py(ro.r("results"))

    #with localconverter(ro.default_converter + pandas2ri.converter):
    #    res_df = ro.conversion.rpy2py(top)

    #####################
    # Ensure a 'gene' column exists
    #if res_df.index.name is None:
    #    res_df.index.name = 'gene'
    #res_df = res_df.reset_index().rename(columns={'index': 'gene'})

    #res_df = res_df.rename_axis('gene').reset_index()

    # Fix gene column
    #if "gene" in res_df.columns:
    #    res_df.rename(columns={"Gene": "gene"}, inplace=True)
    #else:
    #    res_df.index.name = "gene"
    #    res_df = res_df.reset_index()

    # 1. If limma stored the genelist, it will be in the column "ID"
    if "ID" in res_df.columns:
        res_df.rename(columns={"ID": "gene"}, inplace=True)

    # 2. Ensure gene is first column
    gene_col = res_df.pop("gene")
    res_df.insert(0, "gene", gene_col)

    # 3. Remove numeric row index
    res_df.reset_index(drop=True, inplace=True)

    # Save results
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    res_df.to_csv(output_path, index=False)
    print(f"[limma] Saved results to {output_path}")

    # ---------------- Diagnostic plots in Python ----------------
    if plot_prefix is not None:
        # MA plot: x = AveExpr, y = logFC
        if 'AveExpr' in res_df.columns and 'logFC' in res_df.columns:
            plt.figure(figsize=(6,6))
            plt.scatter(res_df['AveExpr'], res_df['logFC'], s=5)
            plt.axhline(0, linestyle='--')
            plt.xlabel('Average expression (log2)')
            plt.ylabel('logFC')
            plt.title('MA plot')
            ma_path = f"{plot_prefix}_MA.png"
            plt.tight_layout()
            plt.savefig(ma_path, dpi=150)
            plt.close()
            print(f"[limma] Saved MA plot to {ma_path}")
        else:
            print('[limma] Cannot make MA plot: AveExpr/logFC not in results')


        # Mean-vs-variance trend: compute gene-wise variance on log_expr
        gene_means = log_expr.mean(axis=1)
        gene_vars = log_expr.var(axis=1)
        plt.figure(figsize=(6,6))
        plt.scatter(gene_means, gene_vars, s=5)
        plt.xlabel('Mean log2 expression')
        plt.ylabel('Variance')
        plt.title('Mean-variance trend (log2 data)')
        mv_path = f"{plot_prefix}_meanvar.png"
        plt.tight_layout()
        plt.savefig(mv_path, dpi=150)
        plt.close()
        print(f"[limma] Saved mean-variance plot to {mv_path}")



def run_edgeR_pseudobulk(counts_df, metadata_df, gene_names, output_path):
    """
    Run edgeR exactTest DEA using pseudobulk counts.
    counts_df: DataFrame (genes x samples)
    metadata_df: DataFrame with columns ['sample', 'condition']
    """

    print("=== run_edgeR_pseudobulk arguments ===")
    print("\ncounts_df head:")
    print(counts_df.head())
    
    print("\nmetadata_df head:")
    print(metadata_df.head())
    
    print("\ngene_names head:")
    print(gene_names[:5])  # assuming gene_names is a list
    
    print("\noutput_path:")
    print(output_path)
    
    # You can optionally also check for NAs
    print("\nCounts NA summary:", counts_df.isna().sum().sum())
    print("Metadata NA summary:", metadata_df.isna().sum().sum())

    pandas2ri.activate()
    edgeR = importr("edgeR")

    # Fill with 0 the NAs
    counts_df = counts_df.fillna(0).astype(int)

    # Wrap gene names for R
    gene_names_df = pd.DataFrame({"gene_id": gene_names})

    with localconverter(ro.default_converter + pandas2ri.converter):
        r_counts = ro.conversion.py2rpy(counts_df)
        r_meta = ro.conversion.py2rpy(metadata_df)
        r_genes = ro.conversion.py2rpy(gene_names_df)

    # R factor for condition
    r_condition = ro.r("as.factor")(r_meta.rx2("condition"))

    # DGEList
    dge = edgeR.DGEList(counts=r_counts, group=r_condition, genes=r_genes)
    keep = edgeR.filterByExpr(dge, group=r_condition)

    keep_py = ro.conversion.rpy2py(keep)
    if not any(keep_py):
        print("️ Skipping: no genes passed filterByExpr.")
        return

    dge_filtered = edgeR.DGEList(
        counts=r_counts.rx(keep, True),
        group=r_condition,
        genes=r_genes.rx(keep, True),
    )
    dge_norm = edgeR.calcNormFactors(dge_filtered)
    dge_disp = edgeR.estimateDisp(dge_norm)
    de = edgeR.exactTest(dge_disp)
    top = edgeR.topTags(de, n=ro.r("nrow")(dge_disp)).rx2("table")

    with localconverter(ro.default_converter + pandas2ri.converter):
        res_df = ro.conversion.rpy2py(top)

    if "genes" in res_df.columns:
        res_df = res_df.rename(columns={"genes": "gene"})

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    res_df.to_csv(output_path, index=False)
    print(f" Saved DE results to {output_path}")


# ---------------- Main ---------------- #

def main():
    args = parse_args()
    run_name = args.base_dir
    samples = args.samples
    out_base = args.out_dir
    cside_dir = args.cside_dir #rctd gene_expr
    celltypes = args.celltypes
    pseudocount = args.pseudocount
    conditions_map = parse_conditions_map(args.conditions_map)
    condition_order = args.condition_order
    all_area_ct = {}  # {area: {celltype: [sample_dfs]}}

    # Iterate samples
    for sample in samples:
        print(f"\n--- Processing sample {sample} ---")
        if sample not in conditions_map:
            raise ValueError(f"Sample {sample} has no condition in conditions_map!")

        if args.delineation_dir != "NO_DELINEATION":
            delim_path = os.path.join(args.delineation_dir, f"{sample}_manual_delineation.csv")
            df_delineation = pd.read_csv(delim_path, index_col=0)
            if df_delineation.columns[0] != "manual_delineation":
                df_delineation.columns = ["manual_delineation"]
            df_delineation.index.name = "barcode"
            df_delineation.index = df_delineation.index.astype(str).str.strip()
        else:
            df_delineation = None

        # Load weights and normalize per spot
        weights_path = os.path.join(run_name, "cell2location_map", sample, f"{sample}_cell2loc_weights.csv")
        if os.path.exists(weights_path):
            df_weights = pd.read_csv(weights_path, index_col=0)
            df_weights.index.name = "barcode"
            df_weights.index = df_weights.index.astype(str).str.strip()
            # normalize weights per barcode to sum to 1
            df_weights = df_weights.div(df_weights.sum(axis=1) + 1e-6, axis=0)
        else:
            print(f"[WARN] No weights file found for {sample}.")
            df_weights = None

        for ct in celltypes:
            expr_path = os.path.join(
                run_name, "cell2location_map", sample, "gene_expr_ct",
                f"{sample}_{ct}_gene_expr.csv"
            )
            if not os.path.exists(expr_path):
                print(f"️ Missing file for {ct} in {sample}, skipping.")
                continue

            df_expr = pd.read_csv(expr_path, index_col=0).T  # rows=barcodes, cols=genes
            df_expr.index.name = "barcode"
            df_expr.index = df_expr.index.astype(str).str.strip()

            # Normalize gene expression by celltype abundance
            if df_weights is not None:
                weight_cols = [c for c in df_weights.columns if c.endswith(ct)]
                if len(weight_cols) > 0:
                    weight_col = weight_cols[0]
                    weights = df_weights.loc[df_expr.index, weight_col].fillna(0)
                    # multiply gene expression by abundance to get per-cell values
                    df_expr = df_expr.mul(weights, axis=0)
                    # zero low-abundance spots
                    df_expr[weights < 0.05] = 0
                else:
                    print(f"[WARN] No abundance column found for {ct} in {sample}.")
            else:
                print(f"[WARN] No weights loaded, skipping normalization for {ct}.")

            # --- END INSERT ---

            if df_delineation is not None:
                df_merged = df_delineation.merge(df_expr, left_index=True, right_index=True, how='inner')
            else:
                df_merged = df_expr.copy()
                df_merged["manual_delineation"] = "WHOLE_SAMPLE"

            if df_merged.shape[0] == 0:
                print(f"[warn] No overlapping barcodes for sample {sample}, ct {ct} (merge empty).")
                continue
            
            for area, area_df in df_merged.groupby("manual_delineation"):
                if area not in all_area_ct:
                    all_area_ct[area] = {}
                if ct not in all_area_ct[area]:
                    all_area_ct[area][ct] = []
                all_area_ct[area][ct].append((sample, area_df))

    # Now process per area × celltype
    for area, ct_dict in all_area_ct.items():
        for ct, sample_dfs in ct_dict.items():
            print(f"\n>>> Processing Area={area}, CellType={ct}")

            filtered_sample_dfs = []

            for s, df in sample_dfs:

                df_numeric = df.select_dtypes(include=[np.number])

                # Filter celltype based on fraction present >0.1 in ≥5% of spots
                n_total_spots = df_numeric.shape[0]
                print(f"total spots in sample {s} and area {area}:\n{n_total_spots}")

                n_spots_present = (df_numeric > 0.1).any(axis=1).sum()
                print(f"of which celltype {ct} has sufficient presence in:\n{n_spots_present}")

                fraction_present = n_spots_present / n_total_spots
                if fraction_present >= 0.05:
                    filtered_sample_dfs.append((s, df))

            if not filtered_sample_dfs:
                print(f"[filter] Skipping {area}-{ct}: insufficient presence in area.")
                continue

            # 2️⃣ Keep only celltypes present in ALL samples (SKIP THIS)
            #samples_with_ct = [s for s, _ in filtered_sample_dfs]
            #if set(samples_with_ct) != set([s for s, _ in sample_dfs]):
            #    print(f"[filter] Skipping {area}-{ct}: not present in all samples, only in ({samples_with_ct}) with sufficient presence")
            #    continue

            # Concatenate filtered spots
            # --- Per-sample mean expression ---
            expr_list = []
            meta_list = []

            for sample, df in filtered_sample_dfs:
                expr_only = df.drop(columns=['manual_delineation'], errors='ignore')
                expr_only = expr_only.select_dtypes(include=[np.number])

                if expr_only.empty:
                    print(f"[warn] No numeric data for {area}-{ct} in {sample}, skipping.")
                    continue

                mean_expr = expr_only.mean(axis=0)  # mean per gene for this sample

                # 🔹 Filter to genes present in CSIDE reference for this sample-area-celltype
                ct_safe = ct.replace(" ", "_")
                ref_gene_path = os.path.join(
                    cside_dir, sample, area, "gene_expr_ct", ct_safe,
                    f"CSIDE_{ct_safe}_expr.csv"
                )
                if os.path.exists(ref_gene_path):
                    df_ref = pd.read_csv(ref_gene_path, index_col=0)
                    ref_genes = df_ref.index.astype(str)
                    common_genes = mean_expr.index.intersection(ref_genes)
                    if len(common_genes) == 0:
                        print(f"[WARN] No overlapping genes for {sample}-{area}-{ct}, skipping.")
                        continue
                    mean_expr = mean_expr.loc[common_genes]
                else:
                    print(f"[WARN] No CSIDE reference found at {ref_gene_path}, skipping filtering for {sample}-{area}-{ct}.")



                expr_list.append(mean_expr)
                meta_list.append({
                    "sample": sample,
                    "condition": conditions_map[sample]
                })

                # Save each sample's expression table
                out_dir = os.path.join(run_name, "cell2location_map", sample,
                     "gene_expr_ct_mean", area)
                os.makedirs(out_dir, exist_ok=True)
                expr_out_csv = os.path.join(out_dir, f"{sample}_{ct}_gene_expr.csv")

                # Write in desired format: gene,mean_expr
                mean_expr.to_csv(expr_out_csv, header=["x"])
                print(f"Saved per-spot expression for {sample}: {expr_out_csv}")

            if not expr_list:
                print(f"️ Skipping {area}-{ct}: no data")
                continue

            # Keep only genes common to all samples
            common_genes_all = set(expr_list[0].index)
            for e in expr_list[1:]:
                common_genes_all &= set(e.index)
            common_genes_all = sorted(common_genes_all)
            expr_list = [e.loc[common_genes_all] for e in expr_list]

            counts_df = pd.DataFrame(expr_list).T  # genes x samples
            counts_df.columns = [m["sample"] for m in meta_list]
            metadata_df = pd.DataFrame(meta_list)
            # Print top expressed genes (mean across samples)
            print(f"Top genes by mean expression in {area}-{ct}:")
            print(counts_df.mean(axis=1).sort_values(ascending=False).head(10))

            # Skip if <2 conditions
            if metadata_df["condition"].nunique() < 2:
                print(f"️ Skipping {area}-{ct}: only one condition")
                continue

            # --- NEW CHECK: require ≥2 replicates per condition (limma requirement) ---
            cond_counts = metadata_df["condition"].value_counts()

            if cond_counts.min() < 2:
                print(f"[skip] {area}-{ct}: not enough replicates per condition for limma: {cond_counts.to_dict()}")
                continue

            # Ensure numeric and fill NA
            counts_df = counts_df.fillna(0.0).astype(float)

            # --- NEW: Keep only genes expressed in all samples ---
            expressed_mask = (counts_df > 0).all(axis=1)
            counts_df = counts_df.loc[expressed_mask, :]
            if counts_df.shape[0] == 0:
                print(f"[warn] No genes expressed in all samples for {area}-{ct}, skipping")
                continue

            # Output paths
            out_dir = os.path.join(out_base, area)
            os.makedirs(out_dir, exist_ok=True)
            out_csv = os.path.join(out_dir, f"DEA_{area}_{ct}_cell2loc_limma.csv")
            plot_pref = os.path.join(out_dir, f"DEA_{area}_{ct}")


            # Run limma pseudobulk
            try:
                run_limma_pseudobulk(counts_df, metadata_df, out_csv, pseudocount=pseudocount, plot_prefix=plot_pref, condition_order=condition_order)
            except Exception as e:
                print(f"[error] Limma failed for {area}-{ct}: {e}")

            #run_edgeR_pseudobulk(counts_df, metadata_df, counts_df.index.tolist(), out_path)


if __name__ == "__main__":
    main()

