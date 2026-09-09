// modules/deg_area_enrichment.nf
process DEG_AREA_ENRICHMENT {
    tag "GSEA/ORA - per-area"
    label 'process_medium'
    container params.container_enrichment
    publishDir "${params.outdir}/DEA_areas_good", mode: 'copy'

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
