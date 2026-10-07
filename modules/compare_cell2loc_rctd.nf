process COMPARE_CELL2LOC_RCTD {
    tag "compare"
    label 'process_low'
    container params.container_rpy2
    publishDir "${params.outdir}/Plots/Comparisons", mode: 'copy'

    input:
    val  ref_level
    val  samples
    val  conditions_map
    val  condition_order
    path cell2loc_map     // "cell2location_map" dir from CELL2LOC_OWNDATA (+ gene_expr_ct_mean from CT_GENE_EXPR_PERCELLTYPE)
    path rctd_results     // "RCTD" dir from CSIDE
    path cell2loc_degs    // "CT_Gene_expr" dir from CT_GENE_EXPR_PERCELLTYPE
    path cside_degs       // "CSIDE" dir from CSIDE (for DEGs)
    path cell2loc_fe      // "FunctionalEnrichment" dir from CT_GENE_EXPR_ENRICHMENT
    path delineation_dir  // "delineation" directory from EXTRACT_SPATIAL_INPUTS, or []
    val  areas            // params.areas list, or [] for whole-sample

    output:
    path "**"

    script:
    // When areas is empty/null we use WHOLE_SAMPLE sentinel so compare scripts
    // apply their whole-sample fallback (no area subdirectory).
    def area_list = areas ? areas : ["WHOLE_SAMPLE"]
    // scCODA differential abundance is per area, so it needs the delineation CSVs
    def da_cmd = delineation_dir ?
        """
        python3 /usr/local/bin/differential_abundance_.py \\
            --input_cell2loc  ${cell2loc_map} \\
            --input_rctd      ${rctd_results} \\
            --delineation_dir ${delineation_dir} \\
            --samples         ${samples.join(' ')} \\
            --conditions_map  ${conditions_map.join(' ')} \\
            --condition_order ${condition_order.join(' ')} \\
            --out_dir         DA
        """ : "echo 'No delineation configured: skipping differential abundance'"
    """
    ${da_cmd}

    python3 /usr/local/bin/compare_cell2loc_rctd.py \\
        --ref_level       ${ref_level} \\
        --output_base_dir . \\
        --pie_chart_dir   PieCharts \\
        --samples         ${samples.join(' ')}

    for AREA in ${area_list.join(' ')}; do
        python3 /usr/local/bin/compare_gene_expr.py \\
            --area         \$AREA \\
            --cell2loc_dir ${cell2loc_map} \\
            --rctd_dir     ${rctd_results} \\
            --output_dir   Gene_Expr/\$AREA

        python3 /usr/local/bin/compare_DEGs.py \\
            --area         \$AREA \\
            --cell2loc_dir ${cell2loc_degs} \\
            --rctd_dir     ${cside_degs} \\
            --output_dir   DEGs/\$AREA

        # NOTE: no enrichment is run on CSIDE DEGs yet, so the CSIDE side is
        # missing and compare_FE.py only warns ("Missing CSIDE FE dir").
        python3 /usr/local/bin/compare_FE.py \\
            --area         \$AREA \\
            --cell2loc_dir ${cell2loc_fe} \\
            --rctd_dir     ${cside_degs}/FunctionalEnrichment \\
            --output_dir   FEs/\$AREA
    done

    python3 /usr/local/bin/plot_global_FE_heatmaps.py \\
        --input_base FEs \\
        --output_dir FEs/global_heatmaps
    """

    stub:
    """
    mkdir -p DA Gene_Expr/WHOLE_SAMPLE DEGs/WHOLE_SAMPLE FEs/WHOLE_SAMPLE FEs/global_heatmaps
    touch DA/placeholder
    """
}
