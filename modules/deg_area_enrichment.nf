// modules/deg_area_enrichment.nf
process DEG_AREA_ENRICHMENT {
    tag "${cid} GSEA/ORA - per-area"
    label 'process_medium'
    container params.container_enrichment
    publishDir "${params.outdir}/contrasts/${cid}/DEA_areas_good", mode: 'copy'

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
    touch FunctionalEnrichment/GSEA/placeholder FunctionalEnrichment/ORA/placeholder
    """
}
