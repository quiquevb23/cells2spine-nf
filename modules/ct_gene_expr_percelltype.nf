process CT_GENE_EXPR_PERCELLTYPE {
    tag "${ref_level}"
    label 'process_medium'
    container params.container_rpy2
    publishDir "${params.outdir}/CT_Gene_expr", mode: 'copy'

    input:
    val  ref_level
    val  masked_celltypes
    val  conditions
    val  samples
    val  manual_celltypes
    val  conditions_map
    val  condition_order
    path cside_results        // CSIDE/** output directory
    path delineation_dir      // "delineation" directory from EXTRACT_SPATIAL_INPUTS, or []

    output:
    path "CT_Gene_expr/**", emit: results

    script:
    def delin_arg = delineation_dir ? "--delineation_dir ${delineation_dir}" : "--delineation_dir NO_DELINEATION"
    """
    python3 /usr/local/bin/ct_gene_expr.py \\
        --base_dir          . \\
        --out_dir           CT_Gene_expr \\
        --cside_dir         ${cside_results} \\
        --ref_level         ${ref_level} \\
        --samples           ${samples.join(' ')} \\
        --celltypes         ${manual_celltypes.join(' ')} \\
        --conditions_map    ${conditions_map.join(' ')} \\
        --condition_order   ${condition_order.join(' ')} \\
        --masked_celltypes  "${masked_celltypes}" \\
        ${delin_arg}
    """

    stub:
    """
    mkdir -p CT_Gene_expr/stub
    touch CT_Gene_expr/stub/ct_gene_expr.csv
    """
}
