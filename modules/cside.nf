process CSIDE {
    tag "${ref_level}"
    label 'process_high'
    container params.container_cside
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val  ref_level
    val  ref_label
    val  masked_celltypes
    val  conditions
    val  samples
    val  conditions_map
    val  condition_order
    val  regions              // params.areas list, or [] for whole-sample
    path cell2loc_map, stageAs: 'cell2location_map/*'  // collected per-sample dirs from CELL2LOC_MAP
    path coords               // "coords" directory from EXTRACT_SPATIAL_INPUTS
    path delineation          // "delineation" directory, or empty when not configured
    path sc_ref_data          // collected reference_<cond> dirs from CELL2LOC_REFERENCE, staged in the work dir: CSIDE.R reads <output_base_dir>/reference_<cond>

    output:
    path "RCTD",  emit: rctd_results
    path "CSIDE", emit: cside_results

    script:
    def regions_arg = regions ? "--regions ${regions.join(',')}" : "--regions ALL_SPOTS"
    def delin_arg   = delineation ? "--delineation_dir ${delineation}" : "--delineation_dir NO_DELINEATION"
    """
    Rscript /usr/local/bin/CSIDE.R \\
        --ref_level         ${ref_level} \\
        --cell2loc_dir      cell2location_map \\
        --coords_dir        ${coords} \\
        ${delin_arg} \\
        --output_base_dir   . \\
        --ref_label         ${ref_label} \\
        --masked_celltypes  "${masked_celltypes}" \\
        --conditions        ${conditions} \\
        --samples           ${samples.join(',')} \\
        --conditions_map    ${conditions_map.join(',')} \\
        --condition_order   ${condition_order.join(',')} \\
        ${regions_arg}
    """

    stub:
    """
    mkdir -p RCTD/stub CSIDE/stub
    touch RCTD/stub/myRCTD.rds CSIDE/stub/cside_result.rds
    """
}
