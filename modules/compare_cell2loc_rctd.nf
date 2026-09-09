// modules/compare_cell2loc_rctd.nf
process COMPARE_CELL2LOC_RCTD {
    tag "compare"
    label 'process_low'
    container params.container_rpy2
    publishDir "${params.outdir}/Plots/Comparisons", mode: 'copy'

    input:
    path cell2loc_map; path rctd_results; path cell2loc_degs; path cside_degs; val areas

    output:
    path "**"

    script:
    """
    differential_abundance_.py --input_cell2loc ${cell2loc_map} --input_rctd ${rctd_results} --out_dir DA
    compare_cell2loc_rctd.py --output_base_dir .
    for AREA in ${areas.join(' ')}; do
        compare_gene_expr.py --area \$AREA --cell2loc_dir ${cell2loc_map} --rctd_dir ${rctd_results} --output_dir Gene_Expr/\$AREA
        compare_DEGs.py --area \$AREA --cell2loc_dir ${cell2loc_degs} --rctd_dir ${cside_degs} --output_dir DEGs/\$AREA
        compare_FE.py --area \$AREA --cell2loc_dir ${cell2loc_degs}/FunctionalEnrichment --rctd_dir ${cside_degs}/FunctionalEnrichment --output_dir FEs/\$AREA
    done
    plot_global_FE_heatmaps.py --input_base FEs --output_dir FEs/global_heatmaps
    """
}
