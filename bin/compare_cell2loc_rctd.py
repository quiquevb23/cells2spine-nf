#!/usr/bin/env python3
'''
This code will do pie charts for each spatial adata to compare the loadings of 
cell2location and RCTD, deconvoluted with the same reference
'''


import os
import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
import seaborn as sns
import matplotlib.patches as mpatches
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Compare and visualize cell2loc vs RCTD loadings.")
    parser.add_argument("--ref_level", type=str, required=True, help="Reference level (e.g., l3)")
# We need to provide: csv_paths of rctd cell2loc path of adata and output path for plots
    parser.add_argument("--output_base_dir", type=str, required=True, help="Base dir containing method outputs")
    parser.add_argument("--pie_chart_dir", type=str, required=True, help="Final output dir for plots, use this")
    #parser.add_argument("--samples", type=str, required=True, help="Sample names")
    parser.add_argument("--samples", type=str, nargs="+", required=True, help="Sample names")

    parser.add_argument("--top_n", type=int, default=20, help="Number of top cell types to plot (default=20)")
    return parser.parse_args()

def plot_pie_overlay(adata, proportions, cell_types, title, output_path, cell_colors=None, border_space=0.05, pie_radius=50, alpha_img=0.5, bw=False):
    # If no cell_colors dict passed, create one internally (fallback)
    if cell_colors is None:
        cell_colors = {cell: sns.color_palette("tab10")[i % 10] for i, cell in enumerate(cell_types)}

    library_id = list(adata.uns["spatial"].keys())[0]

    images = adata.uns["spatial"][library_id]["images"]

    if "hires" in images:
        img_key = "hires"
    else:
        img_key = "lowres"

    img = images[img_key]

    pos = adata.obsm["spatial"].copy()
    pos[:, 1] = -pos[:, 1]

    # Add border space
    x_min, x_max = pos[:, 0].min(), pos[:, 0].max()
    y_min, y_max = pos[:, 1].min(), pos[:, 1].max()
    border_x = (x_max - x_min) * border_space
    border_y = (y_max - y_min) * border_space
    pos[:, 0] = np.where(pos[:, 0] == x_min, x_min - border_x, pos[:, 0])
    pos[:, 0] = np.where(pos[:, 0] == x_max, x_max + border_x, pos[:, 0])
    pos[:, 1] = np.where(pos[:, 1] == y_min, y_min - border_y, pos[:, 1])
    pos[:, 1] = np.where(pos[:, 1] == y_max, y_max + border_y, pos[:, 1])

    # Crop image for visualization
    scalefactors = adata.uns["spatial"][library_id]["scalefactors"]
    
    if img_key == "hires":
        scale_factor = scalefactors["tissue_hires_scalef"]
    else:
        scale_factor = scalefactors["tissue_lowres_scalef"]

    spots = adata.obsm["spatial"] * scale_factor
    x_min, y_min = np.min(spots, axis=0)
    x_max, y_max = np.max(spots, axis=0)
    x_min -= (x_max - x_min) * border_space
    x_max += (x_max - x_min) * border_space
    y_min -= (y_max - y_min) * border_space
    y_max += (y_max - y_min) * border_space

    img_h, img_w = img.shape[:2]
    x_min, x_max = int(np.clip(x_min, 0, img_w - 1)), int(np.clip(x_max, 0, img_w - 1))
    y_min, y_max = int(np.clip(y_min, 0, img_h - 1)), int(np.clip(y_max, 0, img_h - 1))
    cropped_img = img[y_min:y_max, x_min:x_max]

    fig, ax = plt.subplots(figsize=(10, 15))
    ax.imshow(cropped_img, cmap="gray" if bw else None, alpha=alpha_img,
              extent=[pos[:, 0].min(), pos[:, 0].max(), pos[:, 1].min(), pos[:, 1].max()])

    # Define consistent color mapping
    #cell_colors = {cell: sns.color_palette("tab10")[i % 10] for i, cell in enumerate(cell_types)}

    # Draw pie charts
    for (x, y), props in zip(pos, proportions.to_numpy()):
        total = sum(props)
        if total == 0: continue
        start_angle = 0
        for cell_type, proportion in zip(cell_types, props):
            if proportion > 0:
                angle = proportion * 360 / total
                wedge = Wedge((x, y), pie_radius, start_angle, start_angle + angle, color=cell_colors[cell_type])
                ax.add_patch(wedge)
                start_angle += angle

    # Add legend
    legend_patches = [mpatches.Patch(color=cell_colors[cell_type], label=cell_type) for cell_type in cell_types]
    ax.legend(handles=legend_patches, title="Cell Types", loc="center left", bbox_to_anchor=(1.05, 0.5))
    ax.set_title(title)
    ax.set_xticks([]), ax.set_yticks([])
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Saved plot: {output_path}")

def plot_stacked_bar_all_samples(area, method_name, sample_data, top_n, matching_columns, cell_colors, output_dir):
    """
    Plots one stacked bar per sample, grouped by area and method.
    sample_data: dict of sample_name -> Series of average proportions
    """
    samples = sorted(sample_data.keys(), key=lambda s: int(''.join(filter(str.isdigit, s))))
    n_samples = len(samples)
    bar_width = 0.6
    ind = np.arange(n_samples)
    '''
    fig, ax = plt.subplots(figsize=(10, 6))

    bottoms = np.zeros(n_samples)

    for cell_type in matching_columns:
        values = [sample_data[sample].get(cell_type, 0) for sample in samples]
        ax.bar(ind, values, bottom=bottoms, color=cell_colors[cell_type], label=cell_type)
        bottoms += values
    '''
    # Compute average across all samples to determine top N
    mean_props = pd.DataFrame(sample_data).T[matching_columns].mean()
    top_cell_types = mean_props.nlargest(top_n).index.tolist()

    # Remaining go to 'Other'
    grouped_columns = top_cell_types + ['Other']

    # Generate colors
    palette = sns.color_palette("tab20", n_colors=len(top_cell_types))
    updated_cell_colors = {cell: palette[i] for i, cell in enumerate(top_cell_types)}
    updated_cell_colors['Other'] = (0.7, 0.7, 0.7)  # Grey for "Other"

    fig, ax = plt.subplots(figsize=(10, 6))

    bottoms = np.zeros(n_samples)

    for cell_type in grouped_columns:
        values = []
        for sample in samples:
            avg_props = sample_data[sample]
            if cell_type == 'Other':
                val = avg_props.drop(top_cell_types, errors='ignore').sum()
            else:
                val = avg_props.get(cell_type, 0)
            values.append(val)

        ax.bar(ind, values, bottom=bottoms, color=updated_cell_colors[cell_type], label=cell_type)
        bottoms += values

    ax.set_xticks(ind)
    ax.set_xticklabels(samples, rotation=45, ha='right')
    ax.set_ylabel('Average Proportion')
    ax.set_title(f"{area} - {method_name} Stacked Bar")
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
    plt.tight_layout()

    filename = f"{area}_{method_name}_stacked_bar.png"
    output_path = os.path.join(output_dir, filename)
    plt.savefig(output_path, dpi=300)
    plt.close()


def save_method_area_csvs(method_data, matching_columns, output_dir):
    for method_name, area_dict in method_data.items():
        for area, sample_dict in area_dict.items():
            df = pd.DataFrame(index=matching_columns)

            for sample, avg_props in sample_dict.items():
                df[sample] = avg_props.reindex(matching_columns)

            filename = f"{method_name}_{area}_proportions.csv".replace(" ", "_")
            csv_path = os.path.join(output_dir, filename)
            df.to_csv(csv_path)
            print(f"Saved: {csv_path}")

def main():
    args = parse_args()

#Output path is --pie_chart_dir
#Csv path is /storage/gge/Quique/Cells2SpineData/Pilot/spatial/matrices/Deconvolution_owndata/{run_name}/RCTD_weights
#Adatas path is /storage/gge/Quique/Cells2SpineData/Pilot/spatial/matrices/Deconvolution_owndata/{run_name}/cell2location_map

    #cell2loc_adatas_path = os.path.join(args.output_base_dir, "cell2location_map", "cleaned") #here we have the sp{sample}.h5ad, cleaned folder
    cell2loc_adatas_path = os.path.join(args.output_base_dir, "cell2location_map")
    rctd_csv_path = os.path.join(args.output_base_dir, "RCTD") #here we have the folders for each sample with the weights as csv
    output_dir = args.pie_chart_dir
    samples_list = args.samples

    # method_name -> area -> sample -> avg cell proportions
    method_data = {
        "RCTD": defaultdict(dict),
        "cell2loc": defaultdict(dict),
    }

    matching_columns = None  # Set once
    cell_colors = None       # Set once

    for sample in samples_list:
        rctd_path = os.path.join(rctd_csv_path, sample, f"{sample}_RCTD_weights.csv")
        adata_path = os.path.join(cell2loc_adatas_path, f"sp{sample}.h5ad")

        if not os.path.exists(rctd_path):
            print(f"RCTD CSV not found for {sample} at {rctd_path}")
            continue
        if not os.path.exists(adata_path):
            print(f"Adata file not found for {sample} at {adata_path}")
            continue

        # Load data
        rctd = pd.read_csv(rctd_path, index_col=0, encoding="utf-8")
        adata = sc.read_h5ad(adata_path)

        # Match cell types between RCTD and adata
        
        # Cell2loc scores are not loaded in adata.obs, but in .obsm, so extract from there
        #cell2loc_raw_df = adata.obsm['q05_cell_abundance_w_sf']
        cell2loc_raw_df = adata.obsm['means_cell_abundance_w_sf']
        cell2loc_raw = cell2loc_raw_df.to_numpy()
        cell2loc_celltypes = list(adata.uns['mod']['factor_names'])
        
        # Get only the matching celltypes between both
        matching_columns = [ct for ct in rctd.columns if ct in cell2loc_celltypes]
        if not matching_columns:
            print(f"No matching cell types found for sample {sample}")
            continue  # or handle accordingly

        # Create a palette for celltypes
        palette = sns.color_palette("tab20", n_colors=len(matching_columns))
        cell_colors = {cell: palette[i] for i, cell in enumerate(matching_columns)}

        # Create a mapping from celltype name to column index
        celltype_to_index = {ct: i for i, ct in enumerate(cell2loc_celltypes)}

        selected_indices = [celltype_to_index[ct] for ct in matching_columns]

        rctd_l2 = rctd[matching_columns]
        cell2loc_l2 = pd.DataFrame(
            cell2loc_raw[:, selected_indices],
            index=adata.obs_names,
            columns=matching_columns
        )
        
        # Align indices
        shared_barcodes = rctd_l2.index.intersection(cell2loc_l2.index)
        rctd_l2 = rctd_l2.loc[shared_barcodes]
        cell2loc_l2 = cell2loc_l2.loc[shared_barcodes]

        if len(rctd_l2) != len(cell2loc_l2):
            print(f"Mismatch in spot count for {sample}: RCTD={len(rctd_l2)}, Cell2loc={len(cell2loc_l2)}")
            continue
        
        # Normalize
        rctd_norm = rctd_l2.div(rctd_l2.sum(axis=1), axis=0)
        cell2loc_norm = cell2loc_l2.div(cell2loc_l2.sum(axis=1), axis=0)

        # Copy adata and store new loadings
        adata_scores = adata[shared_barcodes].copy()
        for celltype in rctd_norm.columns:
            adata_scores.obs[celltype + '_rctd'] = rctd_norm[celltype].values
        for celltype in cell2loc_norm.columns:
            adata_scores.obs[celltype + '_cell2loc'] = cell2loc_norm[celltype].values

        # Add raw RCTD scores to adata_scores.obs
        rctd_raw = rctd.loc[shared_barcodes, matching_columns]
        for celltype in matching_columns:
            adata_scores.obs[celltype + '_rctd_raw'] = rctd_raw[celltype].values

        # Create directory for raw RCTD spatial plots
        rctd_output_dir = os.path.join(output_dir, "RCTD_deconv_scores")
        os.makedirs(rctd_output_dir, exist_ok=True)

        # Plot spatial plots for raw RCTD scores
        for celltype in matching_columns:
            colname = celltype + '_rctd_raw'
            if colname not in adata_scores.obs.columns:
                continue

            sc.pl.spatial(
                adata_scores,
                cmap='magma',
                color=colname,
                size=1.3,
                img_key=img_key,
                vmin=0,
                vmax='p99.2',
                show=False
            )
            plot_path = os.path.join(rctd_output_dir, f"spatial_RCTD_raw_{sample}_{celltype}.png")
            plt.savefig(plot_path, dpi=300)
            plt.close()
            print(f"Saved spatial plot (raw): {plot_path}")

        # Determine matching columns once (assumes consistent across samples)
        if matching_columns is None:
            matching_columns = rctd_norm.columns.tolist()
            palette = sns.color_palette("tab20", n_colors=len(matching_columns))
            cell_colors = {cell: palette[i] for i, cell in enumerate(matching_columns)}

        # Compute average proportions per area
        for area in adata_scores.obs['manual_delineation'].unique():
            mask = adata_scores.obs['manual_delineation'] == area
            avg_rctd = rctd_norm[mask].mean().reindex(matching_columns)
            avg_cell2loc = cell2loc_norm[mask].mean().reindex(matching_columns)

            method_data["RCTD"][area][sample] = avg_rctd
            method_data["cell2loc"][area][sample] = avg_cell2loc


        # Plot pie chart for RCTD
        plot_pie_overlay(
            adata=adata_scores,
            proportions=rctd_norm,
            cell_types=matching_columns,
            cell_colors=cell_colors,
            title=f"{sample} - RCTD Proportions",
            output_path=os.path.join(output_dir, f"{sample}_piechart_RCTD.png")
        )

        # Plot pie chart for Cell2loc
        plot_pie_overlay(
            adata=adata_scores,
            proportions=cell2loc_norm,
            cell_types=matching_columns,
            cell_colors=cell_colors,
            title=f"{sample} - Cell2location Proportions",
            output_path=os.path.join(output_dir, f"{sample}_piechart_cell2loc.png")
        )
        
    for method_name, method_dict in method_data.items():
        for area, sample_dict in method_dict.items():
            plot_stacked_bar_all_samples(
                area=area,
                method_name=method_name,
                sample_data=sample_dict,
                top_n=args.top_n,
                matching_columns=matching_columns,
                cell_colors=cell_colors,
                output_dir=output_dir
            )

    save_method_area_csvs(method_data, matching_columns, output_dir)

if __name__ == "__main__":
    main()




