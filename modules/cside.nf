// One task per pairwise contrast (see contrasts in main.nf): CSIDE.R and its
// DE design handle exactly two conditions, condition_order = [test, ref].
process CSIDE {
    tag "${cid}"
    label 'process_high'
    container params.container_cside
    publishDir "${params.outdir}/contrasts/${cid}", mode: 'copy'

    input:
    tuple val(cid), val(samples), val(conditions_map), val(condition_order)
    val  ref_level
    val  ref_label
    val  masked_celltypes
    val  regions              // params.areas list, or [] for whole-sample
    path cell2loc_map, stageAs: 'cell2location_map/*'  // collected per-sample dirs from CELL2LOC_MAP
    path coords               // "coords" directory from EXTRACT_SPATIAL_INPUTS
    path delineation          // "delineation" directory, or empty when not configured
    path sc_ref_data          // collected reference_<cond> dirs from CELL2LOC_REFERENCE, staged in the work dir: CSIDE.R reads <output_base_dir>/reference_<cond>

    output:
    tuple val(cid), path("RCTD"),  emit: rctd_results
    tuple val(cid), path("CSIDE"), emit: cside_results
    path "RCTD/*/*_RCTD_weights.csv", emit: rctd_weights, optional: true  // per sample, for maSigPro abundance

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
        --conditions        ${condition_order.join(',')} \\
        --samples           ${samples.join(',')} \\
        --conditions_map    ${conditions_map.join(',')} \\
        --condition_order   ${condition_order.join(',')} \\
        ${regions_arg}
    """

    stub:
    """
    mkdir -p CSIDE/stub
    touch CSIDE/stub/cside_result.rds
    for s in ${samples.join(' ')}; do mkdir -p RCTD/\$s; touch RCTD/\$s/\${s}_RCTD_weights.csv; done
    """
}
