// modules/ct_gene_expr_enrichment.nf
process CT_GENE_EXPR_ENRICHMENT {
    tag "${cid} GSEA/ORA - per-celltype"
    label 'process_medium'
    container params.container_enrichment
    publishDir "${params.outdir}/contrasts/${cid}/CT_Gene_expr", mode: 'copy'

    input:
    tuple val(cid), path(deg_dir)

    output:
    tuple val(cid), path("FunctionalEnrichment"), emit: results

    script:
    """
    mkdir -p FunctionalEnrichment/GSEA FunctionalEnrichment/ORA
    Rscript /usr/local/bin/GSEA.R --deg_dir ${deg_dir} --output_dir FunctionalEnrichment/GSEA
    Rscript /usr/local/bin/ORA.R  --deg_dir ${deg_dir} --output_dir FunctionalEnrichment/ORA
    """

    stub:
    """
    mkdir -p FunctionalEnrichment/GSEA FunctionalEnrichment/ORA
    touch FunctionalEnrichment/GSEA/gsea_results.csv FunctionalEnrichment/ORA/ora_results.csv
    """
}
