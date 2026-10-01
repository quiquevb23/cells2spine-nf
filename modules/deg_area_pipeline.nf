process DEG_AREA_PIPELINE {
    tag "${ref_level}"
    label 'process_medium'
    container params.container_rpy2
    publishDir "${params.outdir}/DEA_areas_good", mode: 'copy'

    input:
    path spatial_input        // raw h5ad root: <spatial_input>/<sample>/outs/matrices/<glob>
    val  h5ad_glob
    val  ref_level
    val  samples
    val  conditions_map
    val  condition_order
    path delineation_dir      // "delineation" directory, or [] when not configured
    val  regions              // params.areas list, or [] for whole-sample

    output:
    path "DEA_areas_good", emit: deg_dir

    script:
    def delin_arg   = delineation_dir ? "--delineation_dir ${delineation_dir}" : "--delineation_dir NO_DELINEATION"
    def regions_arg = regions         ? "--regions ${regions.join(',')}"        : "--regions ALL_SPOTS"
    """
    python3 /usr/local/bin/DEG_areas_pipeline.py \\
        --ref_level         ${ref_level} \\
        --output_base_dir   . \\
        --spatial_input     ${spatial_input} \\
        --spatial_h5ad_glob "${h5ad_glob}" \\
        --samples           ${samples.join(' ')} \\
        --conditions_map    ${conditions_map.join(' ')} \\
        --condition_order   ${condition_order.join(' ')} \\
        ${delin_arg} \\
        ${regions_arg}

    python3 /usr/local/bin/VolcanoPlot_DEG_pipeline.py \\
        --output_base_dir . \\
        --samples         ${samples.join(' ')}
    """

    stub:
    """
    mkdir -p DEA_areas_good/stub
    touch DEA_areas_good/stub/placeholder
    """
}
