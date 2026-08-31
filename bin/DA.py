import argparse
import os
import anndata
import numpy as np
import pandas as pd
import pertpy as pt


def run_da_pertpy(input_dir, manual_celltypes=None, output_dir=None):
    # Load all matching samples
    adatas = []
    for f in os.listdir(input_dir):
        if f.startswith('sp') and f.endswith('.h5ad'):
            adata = sc.read_h5ad(os.path.join(input_dir, f))
            # If not using manual list and have a celltype_key to extract
            if celltype_key and not manual_celltypes:
                if celltype_key not in adata.obsm:
                    raise KeyError(f"{celltype_key} not found in obsm of {f}")
                adata.obs['celltype'] = adata.obsm[celltype_key]
            print(f"Loaded {f}, conditions: {adata.obs['condition'].unique()}")
            adatas.append(adata)

    # Get all unique areas
    areas = set()
    for adata in adatas:
        areas.update(adata.obs['manual_delineation'].unique())

    results = {}
    for area in areas:
        print(f"Processing area: {area}")
        # Subset samples for this area
        area_adatas = []
        for adata in adatas:
            subset = adata[adata.obs['manual_delineation'].astype(str) == str(area)].copy()
            # Use raw counts if you want, but for DA on celltype abundance scores, we use obsm below
            subset.X = subset.layers['counts'].copy() if 'counts' in subset.layers else subset.X
            area_adatas.append(subset)

        combined = anndata.concat(area_adatas)
        combined.obs_names_make_unique()
        print(f"Conditions present: {combined.obs['condition'].unique()}")

        # Determine celltypes for this area
        if manual_celltypes:
            celltypes = manual_celltypes
        else:
            # fallback to keys in obsm for celltype abundance scores
            # or use 'celltype' obs column unique values
            obsm_keys = combined.obsm.keys()
            celltypes = list(obsm_keys)
            if 'celltype' in combined.obs:
                celltypes = combined.obs['celltype'].astype(str).unique()

        # Filter celltypes to those present in all adatas' obsm
        filtered_celltypes = []
        for ct in celltypes:
            if all(ct in a.obsm for a in area_adatas):
                filtered_celltypes.append(ct)
            else:
                print(f"Skipping {ct}, not in all samples' obsm")

        if len(filtered_celltypes) == 0:
            print(f"No common celltypes found for area {area}. Skipping.")
            continue

        # Pseudobulk per sample: average celltype abundance scores for this area
        pseudo_data = []
        for ad in area_adatas:
            sample = ad.uns.get('sample_id', ad.uns.get('name', 'unknown_sample'))
            condition = ad.obs['condition'].unique()
            if len(condition) != 1:
                print(f"Warning: Multiple conditions in sample {sample} for area {area}. Skipping sample.")
                continue
            condition = condition[0]

            # Average celltype abundance for the area in this sample
            mean_scores = {}
            for ct in filtered_celltypes:
                # mean over all spots for this celltype score
                mean_scores[ct] = ad.obsm[ct].mean(axis=0) if ad.obsm[ct].ndim > 1 else ad.obsm[ct].mean()
            mean_scores['sample'] = sample
            mean_scores['condition'] = condition
            pseudo_data.append(mean_scores)

        pseudo_df = pd.DataFrame(pseudo_data)

        # Prepare input for pertpy: abundance matrix and group labels
        X = pseudo_df[filtered_celltypes].values
        groups = pseudo_df['condition'].values

        data = pt.data.from_pandas(X, groups)
        test = pt.test.diff_abundance(data, test='permutation', n_perm=1000)

        res_df = pd.DataFrame({
            'celltype': filtered_celltypes,
            'test_statistic': test.stat,
            'pvalue': test.pvalue
        }).sort_values('pvalue')

        print(f"DA results for area {area}:")
        print(res_df)
        results[area] = res_df

        # Optionally save
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            res_df.to_csv(os.path.join(output_dir, f"DA_results_area_{area}.csv"), index=False)

    return results


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--sp_deconv_data_dir', required=True)
    parser.add_argument('--manual_celltypes', nargs='*', default=[])
    parser.add_argument('--da_output_dir')
    args = parser.parse_args()

    results = run_da_pertpy(
        input_dir=args.sp_deconv_data_dir,
        manual_celltypes=args.manual_celltypes if len(args.manual_celltypes) > 0 else None,
        output_dir=args.da_output_dir
    )

if __name__ == "__main__":
    main()

