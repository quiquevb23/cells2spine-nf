#!/usr/bin/env python3
"""Extract a manual-delineation obs column from spatial AnnData into a
per-sample CSV (barcode, cluster_label), matching the format maSigPro's
R script expects. Skipped entirely if no delineation_obs_key configured."""
import argparse
import scanpy as sc
import pandas as pd
from pathlib import Path

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--spatial_input", required=True)
    p.add_argument("--samples", nargs="+", required=True)
    p.add_argument("--obs_key", required=True,
                   help="obs column holding the manual delineation labels")
    p.add_argument("--output_dir", required=True)
    p.add_argument("--spatial_h5ad_glob", default="*_manual_delineation.h5ad",
                   help="Glob pattern for h5ad files within each sample directory")

    args = p.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    for sample in args.samples:
        sample_path = Path(args.spatial_input) / sample / "outs" / "matrices"
        matches = list(sample_path.glob(args.spatial_h5ad_glob))
        if not matches:
            print(f"[WARN] No h5ad matching '{args.spatial_h5ad_glob}' in {sample_path}, skipping {sample}")
            continue
        h5ad_path = matches[0]

        adata = sc.read_h5ad(h5ad_path)
        if args.obs_key not in adata.obs.columns:
            raise KeyError(f"obs key '{args.obs_key}' not found in {sample}.h5ad "
                            f"(available: {list(adata.obs.columns)})")

        df = pd.DataFrame({
            "X": adata.obs_names,
            args.obs_key: adata.obs[args.obs_key].values
        }).rename(columns={args.obs_key: "manual_delineation"})

        df.to_csv(outdir / f"{sample}_manual_delineation.csv", index=False)
        print(f"[OK] wrote {sample}_manual_delineation.csv ({len(df)} spots)")

if __name__ == "__main__":
    main()
