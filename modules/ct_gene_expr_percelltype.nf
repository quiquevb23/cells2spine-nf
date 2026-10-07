process CT_GENE_EXPR_PERCELLTYPE {
    tag "${ref_level}"
    label 'process_medium'
    container params.container_rpy2
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val  ref_level
    val  masked_celltypes
    val  conditions
    val  samples
    val  manual_celltypes
    val  conditions_map
    val  condition_order
    path cell2loc_map         // "cell2location_map" dir; ct_gene_expr.py reads ./cell2location_map/<sample>/ (base_dir .)
    path cside_results        // "CSIDE" output directory
    path delineation_dir      // "delineation" directory from EXTRACT_SPATIAL_INPUTS, or []

    output:
    path "CT_Gene_expr", emit: results

    script:
    def delin_arg = delineation_dir ? "--delineation_dir ${delineation_dir}" : "--delineation_dir NO_DELINEATION"
    // ct_gene_expr.py only creates CT_Gene_expr/<area> when limma runs (needs >=2 replicates
    // per condition), so make sure the output dir exists even when everything is skipped.
    // NOTE: it also writes cell2location_map/<sample>/gene_expr_ct_mean/ through the staged
    // symlink, i.e. into CELL2LOC_OWNDATA's work dir; COMPARE_CELL2LOC_RCTD reads it from there.
    """
    mkdir -p CT_Gene_expr
    python3 /usr/local/bin/ct_gene_expr.py \\
        --base_dir          . \\
        --out_dir           CT_Gene_expr \\
        --cside_dir         ${cside_results} \\
        --samples           ${samples.join(' ')} \\
        --celltypes         ${manual_celltypes.join(' ')} \\
        --conditions_map    ${conditions_map.join(' ')} \\
        --condition_order   ${condition_order.join(' ')} \\
        ${delin_arg}
    """

    stub:
    """
    mkdir -p CT_Gene_expr/stub
    touch CT_Gene_expr/stub/ct_gene_expr.csv
    """
}
