process COMPARE_CELL2LOC_RCTD {
    tag "${cid}"
    // scCODA HMC (20k samples per area x method) keeps per-spot traces in RAM:
    // OOM-killed at 8 GB, so medium (8 CPU / 32 GB / 8h).
    label 'process_medium'
    container params.container_rpy2
    publishDir "${params.outdir}/contrasts/${cid}/Plots/Comparisons", mode: 'copy'

    input:
    // one task per pairwise contrast; all dirs are that contrast's outputs
    tuple val(cid), val(samples), val(conditions_map), val(condition_order),
          path(rctd_results),   // "RCTD" dir from CSIDE
          path(cside_degs),     // "CSIDE" dir from CSIDE (for DEGs)
          path(cell2loc_degs),  // "CT_Gene_expr" dir from CT_GENE_EXPR_PERCELLTYPE
          path(cell2loc_map, stageAs: 'cell2location_map'),  // CT_GENE_EXPR_PERCELLTYPE's cell2location_map (per-sample dirs + gene_expr_ct_mean)
          path(cell2loc_fe)     // "FunctionalEnrichment" dir from CT_GENE_EXPR_ENRICHMENT
    val  ref_level
    path delineation_dir  // "delineation" directory from EXTRACT_SPATIAL_INPUTS, or []
    val  areas            // params.areas list, or [] for whole-sample

    output:
    // Only this step's own outputs: "**" also re-published every staged input
    // (cell2location_map, RCTD, CSIDE, delineation, models...) into Plots/Comparisons.
    // DEGs/FEs are absent when there are no replicates (empty CT_Gene_expr).
    path "DA",        optional: true
    path "PieCharts", optional: true
    path "Gene_Expr", optional: true
    path "DEGs",      optional: true
    path "FEs",       optional: true

    script:
    // When areas is empty/null we use WHOLE_SAMPLE sentinel so compare scripts
    // apply their whole-sample fallback (no area subdirectory).
    def area_list = areas ? areas : ["WHOLE_SAMPLE"]
    // scCODA differential abundance is per area, so it needs the delineation CSVs
    def da_cmd = delineation_dir ?
        """
        python3 /usr/local/bin/differential_abundance_.py \\
            --input_cell2loc  cell2location_map \\
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
            --cell2loc_dir cell2location_map \\
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
