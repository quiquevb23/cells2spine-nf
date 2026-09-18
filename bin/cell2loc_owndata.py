#!/usr/bin/env python3
'''
	Deconvolution of ST pilot samples with our own annotated data
	- First with celltypes_l2, then with celltypes_l3
'''
import os
import pandas as pd
from sklearn import metrics
import multiprocessing as mp
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import scanpy as sc
import numpy as np
import cell2location
from scipy import sparse
from scipy.sparse import csr_matrix
import re

def parse_args():
    parser = argparse.ArgumentParser(description="Process directories for single-cell data.")

    parser.add_argument('--ref_level', type=str, required=True)
    parser.add_argument('--single_cell_ref', type=str, required=True)
    parser.add_argument('--spatial_input', type=str, required=True)  # renamed from base_dir
    parser.add_argument('--output_base_dir', type=str, required=True)
    parser.add_argument("--conditions", help="Comma-separated list of conditions")
    parser.add_argument('--ref_label', type=str, required=True)
    parser.add_argument('--masked_celltypes', type=str, default="")
    parser.add_argument('--cutoff_celltypes', type=float, default=0.5)

    # NEW: samples (multiple values allowed)
    parser.add_argument("--samples", nargs="+", help="List of spatial sample names")
    parser.add_argument("--conditions_map", type=str, nargs="+", required=True,
        help="Mapping of samples to conditions, in the form 'sample:condition'.")

    # Parse the arguments
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


args = parse_args()

spatial_input = args.spatial_input
output_base_dir = args.output_base_dir
ref_level = args.ref_level
ref_label = args.ref_label
conditions = [c.strip() for c in args.conditions.split(",") if c.strip()]
single_cell_ref_h5 = args.single_cell_ref
masked_celltypes_raw = [x.strip() for x in args.masked_celltypes.split(",") if x.strip()]
cutoff_ct = args.cutoff_celltypes
samples = args.samples
conditions_map = parse_conditions_map(args.conditions_map)


run_name = os.path.join(output_base_dir, "cell2location_map")
os.makedirs(run_name, exist_ok=True)

def sanitize_label(label):
    return re.sub(r'[\\/:*?"<>|]', '_', label)


masked_celltypes = [sanitize_label(ct) for ct in masked_celltypes_raw]

def clean_ct_name(ct):
    return re.sub(r'^q05cell_abundance_w_sf_', '', ct)


def save_cell2loc_weights(adata, sample_out_dir, sample_name):
    cell2loc_raw_df = adata.obsm['q05_cell_abundance_w_sf'] #they recommend this
    cell2loc_celltypes = list(adata.uns['mod']['factor_names'])

    # Build DataFrame with spots as index
    cell2loc_df = pd.DataFrame(
        cell2loc_raw_df.to_numpy(),
        index=cell2loc_raw_df.index,
        columns=cell2loc_celltypes
    )

    # Clean cell type names
    cell2loc_df.columns = [clean_ct_name(ct) for ct in cell2loc_df.columns]
    #save weights in ${RUN_BASE_DIR}/cell2location_map/sample/f"{sample}_cell2loc_weights.csv")
    # Output path
    out_file = os.path.join(sample_out_dir, f"{sample_name}_cell2loc_weights.csv")
    cell2loc_df.to_csv(out_file)



def create_reference_model(adata_ref, condition):
    # Save reference models in 'cell2location_map' folder
    ref_signatures = os.path.join(run_name, f"reference_signatures_{condition}")

    if os.path.exists(ref_signatures) and os.path.exists(os.path.join(ref_signatures, "sc.h5ad")):
        print(f"Reference model for {condition} already exists, skipping.")
        return

    os.makedirs(ref_signatures, exist_ok=True)
    adata_ref.X = adata_ref.layers["counts"].copy()
    adata_ref.var_names_make_unique()
    adata_ref.obs_names_make_unique()
    del adata_ref.raw

    # -----------------------------
    # NEW: mitochondrial filtering
    # -----------------------------
    mt_mask = adata_ref.var_names.str.upper().str.startswith("MT-")
    if mt_mask.any():
        adata_ref = adata_ref[:, ~mt_mask].copy()

    from cell2location.utils.filtering import filter_genes
    #We increase it from 3 to 5% to remove more 'general markers'
    selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.15, nonz_mean_cutoff=1.2)
    adata_ref = adata_ref[:, selected].copy()

    cell2location.models.RegressionModel.setup_anndata(adata_ref, batch_key='sample', labels_key=ref_label)
    mod = cell2location.models.RegressionModel(adata_ref)
    mod.view_anndata_setup()
    mod.train(max_epochs=350, accelerator='cpu')

    # Save ELBO training curve
    history = mod.history
    plt.figure(figsize=(6, 4))
    plt.plot(history["elbo_train"], label="Train ELBO", color="blue")
    if "elbo_val" in history:
        plt.plot(history["elbo_val"], label="Validation ELBO", color="orange")
    plt.xlabel("Epoch")
    plt.ylabel("ELBO")
    plt.title(f"ELBO loss curve - {condition}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(ref_signatures, f"ELBO_loss_{condition}.png"), dpi=300)
    plt.close()

    # Option to select full posterior 1 or only quantiles 2
    #adata_ref = mod.export_posterior(adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500})
    adata_ref = mod.export_posterior(
        adata_ref, use_quantiles=True,
        # choose quantiles
        add_to_varm=["q05","q50", "q95", "q0001"], #only quantiles exported
        #    sample_kwargs={'batch_size': 2500} #full posterior
    )

    mod.save(ref_signatures, overwrite=True)
    adata_ref.write(os.path.join(ref_signatures, "sc.h5ad"))

def deconvolution(adata_vis, sample, condition, cutoff_ct):
    # Run deconvolution and save in run_folders:
    # cell2location_map/{sample}/sp{sample}.h5ad: adata after deconv (general interest)

    # cell2location_map/{sample}/{sample}_cell2loc_weights.csv: weights of celltype abundances (for comparison)

    # cell2location_map/{sample}/gene_expr_ct/{sample}_{ct}_gene_expr.csv: csv for each cell type in each area for specific gene expr (used later for testing with edgeR pseudobulk)

    # output_base_dir/Plots/Cell2loc_weights/{sample}: png for weight in spatial (also save the RCTD weights here)

    #THIS 3 are extracted in conversor_reference.py and conversor.py prior to running RCTD + CSIDE and independently of Cell2loc
    # cell2location_map/{sample}/{sample}_counts.csv: raw spatial counts of celltype abundances (for building RCTD object)
    # cell2location_map/{sample}/{sample}_coordinates.csv: spatial coordinates (for building RCTD object)
    # cell2location_map/{sample}/{sample}_manual_delineation.csv: manual_delineation (for building RCTD object)

    # Define sample-output dir for png:
    plots_dir = os.path.join(output_base_dir, 'Plots', 'Cell2loc', sample)
    os.makedirs(plots_dir, exist_ok=True)

    # Define sample-output dir for weights and per-ct-expr
    sample_out_dir = os.path.join(run_name, sample)
    os.makedirs(sample_out_dir, exist_ok=True)

    adata_vis = adata_vis.copy()


    # Select source matrix
    if 'counts' in adata_vis.layers:
        adata_vis.X = adata_vis.layers['counts'].copy()
    else:
        print("WARNING: 'counts' layer not found, using adata_vis.X instead")

    # Use counts if available, otherwise X
    counts = adata_vis.layers.get('counts', adata_vis.X).T

    # Ensure sparse
    counts = csr_matrix(counts)

    # Save to CSV
    df = pd.DataFrame.sparse.from_spmatrix(
        counts,
        index=adata_vis.var_names,
        columns=adata_vis.obs_names
    )
    df.to_csv(os.path.join(sample_out_dir, f"{sample_name}_counts.csv"))

    # Ensure X is sparse
    if not sparse.issparse(adata_vis.X):
        adata_vis.X = sparse.csr_matrix(adata_vis.X)

    # Check if 'mt' column exists; if not, compute QC metrics incl. mt genes
    # ---------------------------
    if 'mt' not in adata_vis.var.columns:
        # Identify mitochondrial genes automatically if human/mouse gene symbols present
        mt_mask = adata_vis.var_names.str.upper().str.startswith('MT-')
        adata_vis.var['mt'] = mt_mask

        # Run Scanpy QC metrics to populate mt gene counts & fractions
        sc.pp.calculate_qc_metrics(
            adata_vis,
            qc_vars=['mt'],     # tells scanpy this is the mitochondrial gene mask
            inplace=True
        )

    mt_genes = adata_vis.var[adata_vis.var['mt'] == True].index
    adata_vis = adata_vis[:, ~adata_vis.var.index.isin(mt_genes)]
    adata_vis.var_names_make_unique()

    ref_signatures_path = os.path.join(run_name, f"reference_signatures_{condition}")
    adata_ref = sc.read_h5ad(os.path.join(ref_signatures_path, "sc.h5ad"))
    mod = cell2location.models.RegressionModel.load(ref_signatures_path, adata_ref)
    adata_ref = mod.export_posterior(adata_ref, use_quantiles=True, add_to_varm=["q05", "q50", "q95", "q0001"])

    # Get average expression signatures per cell type
    if 'means_per_cluster_mu_fg' in adata_ref.varm:
        # Case 1: full posterior was exported
        inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' for i in adata_ref.uns['mod']['factor_names']]].copy()
    elif 'q50_per_cluster_mu_fg' in adata_ref.varm:
        # Case 2: only quantiles available → use posterior median
        inf_aver = adata_ref.varm['q50_per_cluster_mu_fg'][[f'q50_per_cluster_mu_fg_{i}' for i in adata_ref.uns['mod']['factor_names']]].copy()
    else:
        raise ValueError("No valid posterior summaries (means or q50) found in adata_ref.")

    inf_aver.columns = adata_ref.uns['mod']['factor_names']

    '''
    intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()
    '''

    common = adata_vis.var_names.intersection(inf_aver.index)
    adata_vis = adata_vis[:, common].copy()
    inf_aver = inf_aver.loc[common]

    cell2location.models.Cell2location.setup_anndata(adata=adata_vis)

    adata_vis = deconvolution_training(adata_vis, inf_aver, sample, sample_out_dir, plots_dir, cutoff_ct)
    #adata_vis = deconvolution_loading(adata_vis, sample, sample_out_dir)

    #save cell2loc weights
    save_cell2loc_weights(adata_vis, sample_out_dir, sample)

    # Plot the cell type abundances using mean
    adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']

    library_id = list(adata_vis.uns["spatial"].keys())[0]

    images = adata_vis.uns["spatial"][library_id]["images"]

    if "hires" in images:
        img_key = "hires"
    else:
        img_key = "lowres"

    for label in adata_vis.uns['mod']['factor_names']:
        sc.pl.spatial(adata_vis, cmap='magma', color=label, ncols=5, size=1.3, img_key=img_key, vmin=0, vmax='p99.2', show=False, spot_size=100)
        sanitized_label = sanitize_label(label)
        plt.savefig(os.path.join(plots_dir, f"deconvolution_{sample}_{sanitized_label}.png"), dpi=300, bbox_inches='tight') #change output dir
        plt.close()

    from cell2location import run_colocation

    if "sample" not in adata_vis.obs.columns:
        adata_vis.obs["sample"] = sample

    res_dict, adata_vis = run_colocation(
        adata_vis,
        model_name='CoLocatedGroupsSklearnNMF',
        train_args={'n_fact': np.arange(5, 30), 'sample_name_col': 'sample', 'n_restarts': 3},
        model_kwargs={'alpha': 0.01, 'init': 'random', "nmf_kwd_args": {"tol": 0.000001}},
        export_args={'path': f'{sample_out_dir}/CoLocatedComb_{sample}/'}
    )
    res_dict['n_fact12']['mod'].plot_cell_type_loadings()
    plt.savefig(os.path.join(plots_dir, f"heatmap_colocation_fact12_{sample}.png"), dpi=300, bbox_inches='tight') #change output dir

    res_dict['n_fact29']['mod'].plot_cell_type_loadings()
    plt.savefig(os.path.join(plots_dir, f"heatmap_colocation_fact29_{sample}.png"), dpi=300, bbox_inches='tight') #change output dir

    #save dictionary
    for n_fact, content in res_dict.items():
        mod = content['mod']

        for attr_name in dir(mod):
            if attr_name.startswith("_"):
                continue  # skip private/internal attributes

            attr_value = getattr(mod, attr_name)
    
            # Save if it's a numeric matrix (ndarray or DataFrame)
            if isinstance(attr_value, (np.ndarray, pd.DataFrame)):
                # Convert to DataFrame if needed
                if isinstance(attr_value, np.ndarray):
                    df = pd.DataFrame(attr_value)
                else:
                    df = attr_value.copy()

                # Try to use meaningful row/column names if available
                if hasattr(mod, "cell_type_names") and df.shape[1] == len(mod.cell_type_names):
                    df.columns = mod.cell_type_names
                if hasattr(mod, "factor_names") and df.shape[0] == len(mod.factor_names):
                    df.index = mod.factor_names


def deconvolution_training(adata_vis, inf_aver, sample, sample_out_dir, plots_dir, cutoff_ct=0.5):

##################################### TRAIN NEW MODEL 
    mod = cell2location.models.Cell2location(adata_vis, cell_state_df=inf_aver, N_cells_per_location=8, detection_alpha=20)
    mod.view_anndata_setup()
    mod.train(max_epochs=7000, batch_size=None, train_size=1, accelerator='cpu')

    # --- plot ELBO loss ---
    history = mod.history
    plt.figure(figsize=(6,4))
    plt.plot(history['elbo_train'], label='Train ELBO', color='blue')
    if 'elbo_val' in history:
        plt.plot(history['elbo_val'], label='Validation ELBO', color='orange')

    plt.xlabel('Epoch')
    plt.ylabel('ELBO')
    plt.title(f'ELBO loss curve - {sample}')
    plt.legend()
    plt.tight_layout()  # ensures labels are not cut off

    # Save figure
    plt.savefig(os.path.join(plots_dir, f"ELBO_loss_{sample}.png"), dpi=300)
    plt.show()
    plt.close()

    #mod.plot_history(1000)
    #plt.legend(labels=['full data training'])
    #plt.savefig(os.path.join(output_dir, f"ELBO_loss_{sample}.png")) #change outputdir

    # Export posterior
    adata_vis = mod.export_posterior(
        adata_vis,
        sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs},
        add_to_obsm=['means', 'stds', 'q05', 'q95']
    )
    mod.save(f"{sample_out_dir}/model_{sample}", overwrite=True)

    print("\n===== COMPUTING expected expression (q05) =====")

    expected_dict = mod.module.model.compute_expected_per_cell_type(
        #mod.samples["post_sample_means"], mod.adata_manager
        mod.samples["post_sample_q05"], mod.adata_manager
    )

    print("Returned keys:", expected_dict.keys())
    print("expected_dict['mu'] shape:",
      len(expected_dict['mu']),
      "cell types")

################ Zero out ct with means abundance less than 0.5
    celltypes = adata_vis.uns['mod']['factor_names']
    
    print("\n===== SETUP =====")
    print("Number of spots:", adata_vis.n_obs)
    print("Number of genes:", adata_vis.n_vars)
    print("Cell types:", celltypes)
    print("Cutoff (q05):", cutoff_ct)


    # Take means abundances
    abundance =  adata_vis.obsm['q05_cell_abundance_w_sf'].copy()
    abundance.index = adata_vis.obs_names


    print("\n===== DEBUG: RAW abundance summary =====")
    print(abundance.describe())
    
    # Zero out low-confidence cell types per spot

    abundance_filtered = abundance.where(abundance >= cutoff_ct, 0)

    abundance_filtered = abundance_filtered.loc[
        :, (abundance_filtered >= cutoff_ct).any(axis=0)
    ]

    # Store back
    adata_vis.obsm['q05_cell_abundance_w_sf'] = abundance_filtered

    # align columns just in case some were dropped
    #abundance_raw_aligned = abundance[abundance_filtered.columns]

    # boolean mask: True where ct was >0 but now == 0
    zeroed_mask = (abundance > 0) & (abundance_filtered == 0)

    zeroed_celltypes_per_spot = zeroed_mask.apply(
        lambda row: [clean_ct_name(ct) for ct in row.index[row]],
        axis=1
    )
    
    adata_vis.obs['zeroed_celltypes'] = zeroed_celltypes_per_spot.astype(str)
    adata_vis.uns['zeroed_celltypes_per_spot'] = zeroed_celltypes_per_spot.to_dict()
    #zeroed_dict = adata_vis.obs['zeroed_celltypes'].apply(eval).to_dict()
    zeroed_dict = {spot: eval(z) for spot, z in adata_vis.obs['zeroed_celltypes'].items()}

    # Save per cell-type gene expr in cell2location_map/{sample}/gene_expr_ct/{sample}_{ct}_gene_expr.csv
    gene_expr_ct_dir = os.path.join(sample_out_dir, 'gene_expr_ct')
    os.makedirs(gene_expr_ct_dir, exist_ok=True)

    for i, n in enumerate(mod.factor_names_):
        #adata_vis.layers[n] = expected_dict['mu'][i] #this corrupts adata
        #adata_vis.obsm[n] = expected_dict['mu'][i]
        #save the gene expr of cell types
        #counts_ct = adata_vis.obsm[n].T

        mu_ct = expected_dict['mu'][i]  # shape = (n_cells, n_genes)

        counts_ct = pd.DataFrame.sparse.from_spmatrix(
            #counts_ct,
            mu_ct.T,
            index=adata_vis.var_names,
            columns=adata_vis.obs_names
        )

        counts_ct.to_csv(os.path.join(gene_expr_ct_dir, f"{sample}_{n}_gene_expr_unfiltered.csv"))

        # ---- Apply filtering per spot ----
        counts_ct_filt = counts_ct.copy()

        for spot in counts_ct.columns:
            if n in zeroed_dict.get(spot, []):
                counts_ct_filt[spot] = 0  # zero all genes for this spot
        counts_ct_filt.to_csv(os.path.join(gene_expr_ct_dir, f"{sample}_{n}_gene_expr.csv"))
        #df = pd.DataFrame.sparse.from_spmatrix(counts_ct, index=adata_vis.var_names, columns=adata_vis.obs_names)
        #df.to_csv(os.path.join(gene_expr_ct_dir, f"{sample}_{n}_gene_expr.csv"))
        #After filtering add the matrix to adata
        counts_ct_filt_sparse = sparse.csr_matrix(counts_ct_filt.T.values)  
        adata_vis.obsm[n] = counts_ct_filt_sparse

        print(f"{n} stored in adata_vis.obsm (filtered)")

    '''
    abundance_filtered = abundance_filtered.apply(
        lambda row: [clean_ct_name(ct) for ct in row.index[row != 0]],
        axis=1
    )
    '''

    # Add to obs for plotting (tutorial-style)
    # This line crashes
    #adata_vis.obs[adata_vis.uns['mod']['factor_names']] = abundance_filtered
    #adata_vis.obs[abundance_filtered.columns] = abundance_filtered

    # Each cell type becomes a column with numeric abundance values
    adata_vis.obs = pd.concat([adata_vis.obs, abundance_filtered], axis=1)

    print("\n===== DEBUG: masked abundance (per CT sum) =====")
    print(abundance_filtered.sum().sort_values(ascending=False))

    #adata_vis.obsm['means_cell_abundance_w_sf'] = masked_abundance.values

    adata_vis.write(f"{sample_out_dir}/sp{sample}.h5ad")
    return adata_vis


def deconvolution_loading(adata_vis, sample, sample_out_dir):
######################### RE-LOAD MODEL
    
    # ---- LOAD previously trained model instead of creating a new one ----
    model_path = os.path.join(f"{sample_out_dir}", f"model_{sample}")
    if not os.path.exists(model_path):
        raise FileNotFoundError(f"Trained model not found at {model_path}")
    #mod = cell2location.models.Cell2location.load(model_path, 'model.pt')
    mod = cell2location.models.Cell2location.load(model_path, adata_vis)
    #sample already exists, load it
    adata_vis = sc.read_h5ad(os.path.join(sample_out_dir, f"sp{sample}.h5ad"))
    
    return adata_vis


# Function to filter cell types with < 25 cells
def filter_small_celltypes(adata, celltype_col, min_cells=25):
    counts = adata.obs[celltype_col].value_counts()
    valid_celltypes = counts[counts >= min_cells].index
    removed = counts[counts < min_cells]
    if not removed.empty:
        print(f"Removing cell types with fewer than {min_cells} cells: {', '.join(removed.index)}")
    return adata[adata.obs[celltype_col].isin(valid_celltypes)].copy()


# First create the regression mod for reference single cell
##############################
##############################


#COMMENT TO NOT RUN THE MODEL ON SC DATA AGAIN



adata_ref = sc.read(single_cell_ref_h5)

print(conditions)

for cond in conditions:
    print(cond)
    print(ref_label)
    print(masked_celltypes)
    adata_ref_cond = adata_ref[adata_ref.obs['condition'] == cond].copy()
    adata_ref_cond = adata_ref_cond[~adata_ref_cond.obs[ref_label].isin(masked_celltypes)].copy()
    adata_ref_cond = filter_small_celltypes(adata_ref_cond, ref_label, min_cells=10)
    create_reference_model(adata_ref_cond, cond)



########################################################
#Call black replaced to be dynamic

for sample_name in samples:
    sample_path = os.path.join(spatial_input, sample_name, "outs", "matrices")
    #h5ad_file = f"{sample_name}_aligned_crop.h5ad"
    #h5ad_file = f"adata_{sample_name}_orig.h5ad"
    h5ad_file = f"{sample_name}_manual_delineation.h5ad"

    h5ad_path = os.path.join(sample_path, h5ad_file)
    if os.path.exists(h5ad_path):
        adata = sc.read_h5ad(h5ad_path)

        # Use condition from adata.obs or map externally if needed
        #conditions_in_sample = adata.obs['condition'].unique()
        #assert len(conditions_in_sample) == 1, f"Multiple conditions found in {sample_name}"
        #cond = conditions_in_sample[0]
        
        # Map conditions from dictionary
        cond = conditions_map[sample_name]
        deconvolution(adata, sample_name, cond, cutoff_ct)

'''
for file_name in os.listdir(spatial_input):
    if file_name.startswith('Spatial'):
        sample_name = file_name
        sample_path = os.path.join(spatial_input, sample_name, "outs", "matrices")
        h5ad_file = f"adata_{sample_name}_manual_delineation_good.h5ad"
        h5ad_path = os.path.join(sample_path, h5ad_file)
        if os.path.exists(h5ad_path):
            adata = sc.read_h5ad(h5ad_path)
            # Conditions used for pairwise testing: give just 2 groups
            conditions = adata.obs['condition'].unique()
            assert len(conditions) == 1, f"Multiple conditions found in {sample_name}"
            deconvolution(adata, sample_name, conditions[0])
'''
