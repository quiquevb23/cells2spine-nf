#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Module imports
// ---------------------------------------------------------------------------
include { EXTRACT_SPATIAL_INPUTS   } from './modules/extract_spatial_inputs'
include { CELL2LOC_OWNDATA         } from './modules/cell2loc_owndata'
include { CSIDE                    } from './modules/cside'
include { CT_GENE_EXPR_PERCELLTYPE } from './modules/ct_gene_expr_percelltype'
include { CT_GENE_EXPR_ENRICHMENT  } from './modules/ct_gene_expr_enrichment'
include { DEG_AREA_PIPELINE        } from './modules/deg_area_pipeline'
include { DEG_AREA_ENRICHMENT      } from './modules/deg_area_enrichment'
include { MASIGPRO_ANALYSIS        } from './modules/masigpro'
include { COMPARE_CELL2LOC_RCTD    } from './modules/compare_cell2loc_rctd'
// include { VALIDATE_PER_AREA        } from './modules/validate_per_area'      // TODO
// include { VALIDATE_PER_CELLTYPE    } from './modules/validate_per_celltype'  // TODO

// ---------------------------------------------------------------------------
// Parameter defaults (biology-specific values come from --params-file JSON)
// ---------------------------------------------------------------------------
params.outdir             = "./results"
params.ref_level          = null
params.single_cell_ref    = null
params.spatial_input      = null
params.spatial_h5ad_glob  = "*.h5ad"
params.ref_label          = null
params.conditions         = null
params.condition_order    = []
params.cutoff_celltypes   = 0
params.delineation_obs_key = null   // null → no delineation extraction
params.areas              = []      // empty → scripts fall back to WHOLE_SAMPLE / ALL_SPOTS
params.masked_celltypes   = ""
params.manual_celltypes   = []
params.sample_metadata    = []      // [{sample, condition, time}, ...]
params.gene_list          = null    // optional path for maSigPro
params.masigpro           = [:]     // {clusters: [...], gene_lists: [...]}

// ---------------------------------------------------------------------------
// Derived values from sample_metadata
// ---------------------------------------------------------------------------
def samples         = params.sample_metadata.collect { it.sample }
def conditions_map  = params.sample_metadata.collect { "${it.sample}:${it.condition}" }
def timepoints      = params.sample_metadata.collect { it.time }.findAll { it != null }.unique()
def n_timepoints    = timepoints.size()

// maSigPro clusters (optional)
def masigproClusters = params.masigpro?.clusters ?: []

// Area list (empty → whole-sample fallback inside scripts)
def areaList = params.areas ?: []

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------
workflow {

    // ------------------------------------------------------------------
    // 1. Extract spatial inputs (counts, coords, delineation) from h5ad
    // ------------------------------------------------------------------
    EXTRACT_SPATIAL_INPUTS(
        file(params.spatial_input),
        samples,
        params.delineation_obs_key,   // null → no delineation CSV written
        params.spatial_h5ad_glob
    )

    // Channel: delineation directory or empty list when not configured
    def delineationCh = params.delineation_obs_key
        ? EXTRACT_SPATIAL_INPUTS.out.delineation
        : Channel.value([])

    // ------------------------------------------------------------------
    // 2. Cell2Location deconvolution + reference extraction
    // NOTE: cell2loc_owndata.py reads h5ad directly — pass spatial_input,
    //       NOT counts (those are only for downstream CSV-based steps).
    // ------------------------------------------------------------------
    CELL2LOC_OWNDATA(
        params.ref_level,
        file(params.single_cell_ref),
        file(params.spatial_input),
        params.ref_label,
        params.conditions,
        params.masked_celltypes,
        params.cutoff_celltypes,
        params.manual_celltypes,
        samples,
        conditions_map
    )

    // ------------------------------------------------------------------
    // 3. CSIDE (runs RCTD internally, then CSIDE DE per region)
    // ------------------------------------------------------------------
    CSIDE(
        params.ref_level,
        params.ref_label,
        params.masked_celltypes,
        params.conditions,
        samples,
        conditions_map,
        params.condition_order,
        areaList,
        CELL2LOC_OWNDATA.out.cell2loc_map,
        EXTRACT_SPATIAL_INPUTS.out.coords,
        delineationCh
    )

    // ------------------------------------------------------------------
    // 4. Per-celltype gene expression & enrichment
    // ------------------------------------------------------------------
    CT_GENE_EXPR_PERCELLTYPE(
        params.ref_level,
        params.masked_celltypes,
        params.conditions,
        samples,
        params.manual_celltypes,
        conditions_map,
        params.condition_order,
        CSIDE.out.cside_results,
        delineationCh
    )

    CT_GENE_EXPR_ENRICHMENT(
        CT_GENE_EXPR_PERCELLTYPE.out.results
    )

    // ------------------------------------------------------------------
    // 5. DEG by area (independent branch, parallel with CSIDE chain)
    // ------------------------------------------------------------------
    DEG_AREA_PIPELINE(
        EXTRACT_SPATIAL_INPUTS.out.counts,
        params.ref_level,
        samples,
        conditions_map,
        params.condition_order,
        delineationCh,
        areaList
    )

    DEG_AREA_ENRICHMENT(
        DEG_AREA_PIPELINE.out.deg_dir
    )

    // ------------------------------------------------------------------
    // 6. maSigPro (only when ≥3 distinct timepoints)
    // ------------------------------------------------------------------
    if (n_timepoints >= 3) {
        // Write sample_metadata to a JSON/CSV file for the R script
        def metadataFile = file("${params.outdir}/sample_metadata.csv")

        // Build CSV content and stage it
        def metaCsv = "sample,condition,time\n" +
            params.sample_metadata.collect { m ->
                "${m.sample},${m.condition},${m.time}"
            }.join("\n")

        // Write via shell — we need it as a staged path for the process
        def geneListPath = params.gene_list ? file(params.gene_list) : []

        MASIGPRO_ANALYSIS(
            EXTRACT_SPATIAL_INPUTS.out.counts,
            delineationCh,
            Channel.value(writeMetadataCsv(params.sample_metadata)),
            Channel.value(geneListPath),
            masigproClusters
        )
    } else {
        log.info "maSigPro skipped: only ${n_timepoints} distinct timepoint(s) " +
                 "(need ≥ 3). Timepoints found: ${timepoints}"
    }

    // ------------------------------------------------------------------
    // 7. Comparisons (cell2loc vs RCTD/CSIDE)
    // ------------------------------------------------------------------
    COMPARE_CELL2LOC_RCTD(
        CELL2LOC_OWNDATA.out.cell2loc_map,
        CSIDE.out.rctd_results,
        CT_GENE_EXPR_PERCELLTYPE.out.results,
        CSIDE.out.cside_results,
        areaList
    )
}

// ---------------------------------------------------------------------------
// Helper: write sample_metadata list → CSV file staged for a process
// ---------------------------------------------------------------------------
def writeMetadataCsv(List metadata) {
    def f = file("${workDir}/sample_metadata_${workflow.sessionId}.csv")
    def lines = ["sample,condition,time"] +
        metadata.collect { m -> "${m.sample},${m.condition},${m.time}" }
    f.text = lines.join("\n") + "\n"
    return f
}

