process MASIGPRO_ANALYSIS {
    tag "maSigPro"
    label 'process_medium'
    container params.container_enrichment
    publishDir "${params.outdir}/masigpro_results", mode: 'copy'

    input:
    path counts_dir
    path delineation_dir
    path sample_metadata
    path gene_list
    val clusters

    output:
    path "masigpro_results/**", emit: results

    script:
    def cluster_arg = clusters ? "--clusters ${clusters.join(',')}" : "--clusters ALL_SPOTS"
    def gene_arg     = gene_list ? "--gene_list ${gene_list}" : "--gene_list NO_FILE"
    def delin_arg    = delineation_dir ? "--delineation_dir ${delineation_dir}" : "--delineation_dir NO_DELINEATION"
    """
    run_masigpro.R \\
        --counts_dir ${counts_dir} \\
        ${delin_arg} \\
        --sample_metadata ${sample_metadata} \\
        ${gene_arg} \\
        ${cluster_arg} \\
        --outdir masigpro_results
    """
}
