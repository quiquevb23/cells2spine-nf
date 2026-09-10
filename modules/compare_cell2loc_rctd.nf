process COMPARE_CELL2LOC_RCTD {
    tag "compare"
    label 'process_low'
    container params.container_rpy2
    publishDir "${params.outdir}/Plots/Comparisons", mode: 'copy'

    input:
    path cell2loc_map     // CELL2LOC_OWNDATA output
    path rctd_results     // RCTD/** from CSIDE
    path cell2loc_degs    // CT_Gene_expr results from CT_GENE_EXPR_PERCELLTYPE
    path cside_degs       // CSIDE/** from CSIDE (for DEGs)
    val  areas            // params.areas list, or [] for whole-sample

    output:
    path "**"

    script:
    // When areas is empty/null we use WHOLE_SAMPLE sentinel so compare scripts
    // apply their whole-sample fallback (no area subdirectory).
    def area_list = areas ? areas : ["WHOLE_SAMPLE"]
    """
    differential_abundance_.py \\
        --input_cell2loc ${cell2loc_map} \\
        --input_rctd     ${rctd_results} \\
        --out_dir        DA

    compare_cell2loc_rctd.py --output_base_dir .

    for AREA in ${area_list.join(' ')}; do
        compare_gene_expr.py \\
            --area         \$AREA \\
            --cell2loc_dir ${cell2loc_map} \\
            --rctd_dir     ${rctd_results} \\
            --output_dir   Gene_Expr/\$AREA

        compare_DEGs.py \\
            --area         \$AREA \\
            --cell2loc_dir ${cell2loc_degs} \\
            --rctd_dir     ${cside_degs} \\
            --output_dir   DEGs/\$AREA

        compare_FE.py \\
            --area         \$AREA \\
            --cell2loc_dir ${cell2loc_degs}/FunctionalEnrichment \\
            --rctd_dir     ${cside_degs}/FunctionalEnrichment \\
            --output_dir   FEs/\$AREA
    done

    plot_global_FE_heatmaps.py \\
        --input_base FEs \\
        --output_dir FEs/global_heatmaps
    """

    stub:
    """
    mkdir -p DA Gene_Expr/WHOLE_SAMPLE DEGs/WHOLE_SAMPLE FEs/WHOLE_SAMPLE FEs/global_heatmaps
    touch DA/placeholder
    """
}
