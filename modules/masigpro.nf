process MASIGPRO_ANALYSIS {
    tag "maSigPro analysis"
    label 'process_medium'
    container "${params.container_enrichment}"

    input:
    path delineation_dir
    path counts_dir
    path sample_metadata
    path gene_list

    output:
    path "masigpro_results/**", emit: results

    script:
    def gene_arg = gene_list ? "--gene_list ${gene_list}" : "--gene_list NO_FILE"
    """
    run_masigpro.R \\
        --delineation_dir ${delineation_dir} \\
        --counts_dir ${counts_dir} \\
        --sample_metadata ${sample_metadata} \\
        ${gene_arg} \\
        --outdir masigpro_results
    """
}
