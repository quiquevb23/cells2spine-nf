// modules/cell2loc_owndata.nf
//
// cell2location is split so memory does not scale with the number of samples:
//   CELL2LOC_REFERENCE : one job per condition — trains the RegressionModel and
//                        exports data/reference_<cond> (conversor_ref.py) for RCTD/CSIDE
//   CELL2LOC_MAP       : one job per sample — Cell2location mapping with its
//                        condition's reference signatures
// Downstream processes re-assemble the per-sample dirs as ./cell2location_map/<sample>
// and the references as ./data/reference_<cond> via `stageAs`.

process CELL2LOC_REFERENCE {
    tag "${condition}"
    label 'process_high'
    container params.container_cell2loc
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val  condition
    path single_cell_ref
    val  ref_label
    val  masked_celltypes

    output:
    tuple val(condition), path("cell2location_map/reference_signatures_${condition}"), emit: signatures
    path "data/reference_${condition}",                                                emit: sc_ref_data

    script:
    """
    python3 /usr/local/bin/cell2loc_owndata.py \\
        --mode reference \\
        --condition ${condition} \\
        --single_cell_ref ${single_cell_ref} \\
        --ref_label ${ref_label} \\
        --masked_celltypes "${masked_celltypes}" \\
        --output_base_dir .

    python3 /usr/local/bin/conversor_ref.py \\
        --sc_ref_path ${single_cell_ref} --output_base_dir . \\
        --ref_label ${ref_label} --conditions ${condition} \\
        --masked_celltypes "${masked_celltypes}"
    """

    stub:
    """
    mkdir -p cell2location_map/reference_signatures_${condition} data/reference_${condition}
    touch cell2location_map/reference_signatures_${condition}/sc.h5ad
    touch data/reference_${condition}/counts.mtx
    """
}

process CELL2LOC_MAP {
    tag "${sample}"
    label 'process_cell2loc_map'
    container params.container_cell2loc
    publishDir "${params.outdir}", mode: 'copy'

    input:
    // reference dir is staged where cell2loc_owndata.py expects it:
    // ./cell2location_map/reference_signatures_<condition>
    tuple val(condition), val(sample), path(ref_signatures, stageAs: 'cell2location_map/*')
    path spatial_input
    val  h5ad_glob
    val  cutoff_celltypes

    output:
    path "cell2location_map/${sample}", emit: sample_dir   // staged downstream as cell2location_map/<sample>
    path "Plots/Cell2loc/${sample}",    emit: plots   // per-sample path: publishing a shared "Plots" dir would let samples overwrite each other

    script:
    """
    python3 /usr/local/bin/cell2loc_owndata.py \\
        --mode map \\
        --sample ${sample} \\
        --condition ${condition} \\
        --spatial_input ${spatial_input} \\
        --spatial_h5ad_glob "${h5ad_glob}" \\
        --cutoff_celltypes ${cutoff_celltypes} \\
        --output_base_dir .
    """

    stub:
    """
    mkdir -p cell2location_map/${sample} Plots/Cell2loc/${sample}
    touch cell2location_map/${sample}/${sample}_counts.csv
    """
}
