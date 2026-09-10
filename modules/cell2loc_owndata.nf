// modules/cell2loc_owndata.nf
process CELL2LOC_OWNDATA {
    tag "${ref_level}"
    label 'process_high'
    container params.container_cell2loc
    publishDir "${params.outdir}/cell2location_map", mode: 'copy'

    input:
    val ref_level; path single_cell_ref; path spatial_input
    val ref_label; val conditions; val masked_celltypes
    val cutoff_celltypes; val manual_celltypes; val samples; val conditions_map

    output:
    path "cell2location_map/**", emit: cell2loc_map

    script:
    """
    cell2loc_owndata.py \\
        --ref_level ${ref_level} --single_cell_ref ${single_cell_ref} \\
        --spatial_input ${spatial_input} --output_base_dir . \\
        --ref_label ${ref_label} --conditions ${conditions} \\
        --masked_celltypes "${masked_celltypes}" --cutoff_celltypes ${cutoff_celltypes} \\
        --samples ${samples.join(' ')} --conditions_map ${conditions_map.join(' ')}

    conversor_ref.py \\
        --sc_ref_path ${single_cell_ref} --output_base_dir . \\
        --ref_label ${ref_label} --conditions ${conditions} \\
        --masked_celltypes "${masked_celltypes}"
    """

    stub:
    """
    mkdir -p cell2location_map/stub
    touch cell2location_map/stub/cell2loc_counts.csv
    """
}
