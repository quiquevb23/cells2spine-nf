// modules/ct_gene_expr_percelltype.nf
process CT_GENE_EXPR_PERCELLTYPE {
    tag "${ref_level}"
    label 'process_medium'
    container params.container_rpy2

    input:
    val ref_level; val masked_celltypes; val conditions; val samples
    val manual_celltypes; val conditions_map; val condition_order
    path cside_results; path delineation_dir

    output:
    path "CT_Gene_expr/**", emit: deg_dir

    script:
    def delin_arg = delineation_dir ? "--delineation_dir ${delineation_dir}" : "--delineation_dir NO_DELINEATION"
    """
    ct_gene_expr.py \\
        --base_dir . --out_dir CT_Gene_expr --cside_dir . \\
        --samples ${samples.join(' ')} --celltypes ${manual_celltypes.join(' ')} \\
        --conditions_map ${conditions_map.join(' ')} --condition_order ${condition_order.join(' ')} \\
        ${delin_arg}
    """
}
