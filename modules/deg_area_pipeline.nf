// modules/deg_area_pipeline.nf
process DEG_AREA_PIPELINE {
    tag "${ref_level}"
    label 'process_medium'
    container params.container_rpy2
    publishDir "${params.outdir}/DEA_areas_good", mode: 'copy'

    input:
    path counts; path coords; val ref_level; val samples; val conditions_map; val condition_order

    output:
    path "DEA_areas_good/**", emit: deg_dir

    script:
    """
    DEG_areas_pipeline.py \\
        --ref_level ${ref_level} --output_base_dir . \\
        --samples ${samples.join(' ')} --conditions_map ${conditions_map.join(' ')} \\
        --condition_order ${condition_order.join(' ')}

    VolcanoPlot_DEG_pipeline.py --output_base_dir . --samples ${samples.join(' ')}
    """
}
