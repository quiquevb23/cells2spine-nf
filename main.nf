nextflow.enable.dsl = 2

include { EXTRACT_SPATIAL_INPUTS }    from './modules/extract_spatial_inputs'
include { CELL2LOC_OWNDATA }          from './modules/cell2loc_owndata'
include { CSIDE }                     from './modules/cside'
include { CT_GENE_EXPR_PERCELLTYPE }  from './modules/ct_gene_expr_percelltype'
include { CT_GENE_EXPR_ENRICHMENT }   from './modules/ct_gene_expr_enrichment'
include { COMPARE_CELL2LOC_RCTD }     from './modules/compare_cell2loc_rctd'
include { DEG_AREA_PIPELINE }         from './modules/deg_area_pipeline'
include { DEG_AREA_ENRICHMENT }       from './modules/deg_area_enrichment'
include { MASIGPRO_ANALYSIS }         from './modules/masigpro'
// include { VALIDATE_PER_AREA }      from './modules/validate_per_area'
// include { VALIDATE_PER_CELLTYPE }  from './modules/validate_per_celltype'

workflow {

    // ── Derive canonical sample/condition lists from ONE source of truth ──
    def samples         = params.sample_metadata.collect { it.sample }
    def conditions_map  = params.sample_metadata.collect { "${it.sample}:${it.condition}" }
    def n_timepoints    = params.sample_metadata.collect { it.time }.findAll { it != null }.unique().size()
    def delineationCh = params.delineation_obs_key ? EXTRACT_SPATIAL_INPUTS.out.delineation : []
    def areaList       = params.areas ?: []   // empty = let each script fall back to WHOLE_SAMPLE    
    
    // ── Shared upstream extraction: counts + coords (+ optional delineation) ──
    EXTRACT_SPATIAL_INPUTS(
        file(params.spatial_input),
        samples,
        params.delineation_obs_key ?: null
    )

    // ── Main branch: cell2loc → CSIDE → ct_gene_expr → compare ──
    CELL2LOC_OWNDATA(
        params.ref_level, file(params.single_cell_ref), file(params.spatial_input),
        params.ref_label, params.conditions, params.masked_celltypes,
        params.cutoff_celltypes, params.manual_celltypes, samples, conditions_map
    )

    CSIDE(
        params.ref_level, params.ref_label, params.masked_celltypes,
        params.conditions, samples, conditions_map, params.condition_order,
        areaList,
        CELL2LOC_OWNDATA.out.cell2loc_map,
        EXTRACT_SPATIAL_INPUTS.out.coords,
        delineationCh
    )

    CT_GENE_EXPR_PERCELLTYPE(
        params.ref_level, params.masked_celltypes, params.conditions,
        samples, params.manual_celltypes, conditions_map,
        params.condition_order, CSIDE.out.cside_results, delineationCh
    )
    CT_GENE_EXPR_ENRICHMENT(CT_GENE_EXPR_PERCELLTYPE.out.deg_dir)

    COMPARE_CELL2LOC_RCTD(
        CELL2LOC_OWNDATA.out.cell2loc_map,
        CSIDE.out.rctd_results,
        CT_GENE_EXPR_PERCELLTYPE.out.deg_dir,
        CSIDE.out.cside_results,
        params.areas
    )

    // ── Validation: not yet resolved, skipped for now ──
    // VALIDATE_PER_AREA(...)
    // VALIDATE_PER_CELLTYPE(...)

    // ── Independent branch: pairwise per-area DEG ──
    DEG_AREA_PIPELINE(
        params.ref_level, file(params.spatial_input), samples, conditions_map,
        params.condition_order, delineationCh
    )
    DEG_AREA_ENRICHMENT(DEG_AREA_PIPELINE.out.deg_dir)

    // ── maSigPro: independent branch, requires >=3 distinct timepoints ──
    if (n_timepoints >= 3) {
        MASIGPRO_ANALYSIS(
            EXTRACT_SPATIAL_INPUTS.out.counts, delineationCh, metaFile,
            params.masigpro.gene_lists ? file(params.masigpro.gene_lists.join(',')) : [],
            areaList
        )
    } else {
        log.warn "Skipping maSigPro: fewer than 3 distinct timepoints in sample_metadata (found ${n_timepoints})"
    }
}
