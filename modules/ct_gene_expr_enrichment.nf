// modules/ct_gene_expr_enrichment.nf
process CT_GENE_EXPR_ENRICHMENT {
    tag "GSEA/ORA - per-celltype"
    label 'process_medium'
    container params.container_enrichment
    publishDir "${params.outdir}/CT_Gene_expr", mode: 'copy'

    input:
    path deg_dir

    output:
    path "FunctionalEnrichment/**"

    script:
    """
    mkdir -p FunctionalEnrichment/GSEA FunctionalEnrichment/ORA
    Rscript /usr/local/bin/GSEA.R --deg_dir ${deg_dir} --output_dir FunctionalEnrichment/GSEA
    Rscript /usr/local/bin/ORA.R  --deg_dir ${deg_dir} --output_dir FunctionalEnrichment/ORA
    """
}
