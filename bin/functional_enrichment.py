#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import gseapy as gp
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description="Run functional enrichment on DEG lists")
    parser.add_argument("--deg_dir", type=str, required=True,
                        help="Directory containing DEG result files per cell type")
    parser.add_argument("--output_dir", type=str, required=True,
                        help="Directory to save enrichment results")
    parser.add_argument("--manual_celltypes", nargs='*', default=[],
                        help="List of manual cell types to process")
    parser.add_argument("--samples", nargs='*', default=[],
                        help="List of samples to process")
    parser.add_argument("--padj_cutoff", type=float, default=0.05,
                        help="Adjusted p-value cutoff for DEG filtering")
    parser.add_argument("--logfc_cutoff", type=float, default=0.25,
                        help="Absolute log fold change cutoff for DEG filtering")
    return parser.parse_args()

def run_enrichment(gene_list, gene_set, organism="Rat"):
    """
    Run enrichment using GSEApy enrichr.
    gene_set examples: 'GO_Biological_Process_2021', 'KEGG_2021_Mouse'
    """
    try:
        enr = gp.enrichr(gene_list=gene_list,
                         description='enrichment',
                         gene_sets=gene_set,
                         organism=organism,
                         outdir=None,  # We'll save manually
                         cutoff=0.5)
        return enr.results
    except Exception as e:
        print(f"Failed to run enrichment for {gene_set}: {e}")
        return pd.DataFrame()

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    deg_dir = Path(args.deg_dir)

    # Detect DEG files
    deg_files = list(deg_dir.glob("*.csv")) + list(deg_dir.glob("*.tsv"))
    if not deg_files:
        raise FileNotFoundError(f"No DEG files found in {deg_dir}")

    for deg_file in deg_files:
        # Derive identifiers from filename
        fname = deg_file.stem
        matching_celltypes = [ct for ct in args.manual_celltypes if ct in fname]
        matching_samples = [s for s in args.samples if s in fname]
        if args.manual_celltypes and not matching_celltypes:
            continue
        if args.samples and not matching_samples:
            continue

        print(f"Processing DEG file: {deg_file.name}")

        # Read DEGs
        try:
            df = pd.read_csv(deg_file, sep=None, engine="python")
        except Exception as e:
            print(f"Skipping {deg_file.name}, cannot read: {e}")
            continue

        # Expect columns: gene, logFC, p_val_adj
        colnames = {c.lower(): c for c in df.columns}
        required_cols = ["gene", "logfc", "fdr"]
        if not all(rc in colnames for rc in required_cols):
            print(f"Skipping {deg_file.name}, missing required columns.")
            continue

        df_filtered = df[
            (df["FDR"] <= args.padj_cutoff) &
            (df["logFC"].abs() >= args.logfc_cutoff)
        ]
        gene_list = df_filtered["gene"].dropna().astype(str).tolist()


        if not gene_list:
            print(f"No significant DEGs for {deg_file.name}")
            continue

        # Run enrichment for GO BP and KEGG
        for gene_set in ["GO_Biological_Process_2021", "KEGG_2021_Rat"]:
            results = run_enrichment(gene_list, gene_set)
            if results.empty:
                continue

            # Save CSV
            out_csv = Path(args.output_dir) / f"{fname}_{gene_set}_enrichment.csv"
            results.to_csv(out_csv, index=False)

            # Save barplot
            try:
                gp.barplot(results,
                           title=f"{fname} - {gene_set}",
                           ofname=str(Path(args.output_dir) / f"{fname}_{gene_set}_barplot.png"))
            except Exception as e:
                print(f"Failed to plot {gene_set} for {fname}: {e}")

        print(f"Finished enrichment for {fname}")

if __name__ == "__main__":
    main()

