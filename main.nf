#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Module imports
// ---------------------------------------------------------------------------
include { EXTRACT_SPATIAL_INPUTS   } from './modules/extract_spatial_inputs'
include { CELL2LOC_REFERENCE       } from './modules/cell2loc_owndata'
include { CELL2LOC_MAP             } from './modules/cell2loc_owndata'
include { CSIDE                    } from './modules/cside'
include { CT_GENE_EXPR_PERCELLTYPE } from './modules/ct_gene_expr_percelltype'
include { CT_GENE_EXPR_ENRICHMENT  } from './modules/ct_gene_expr_enrichment'
include { DEG_AREA_PIPELINE        } from './modules/deg_area_pipeline'
include { DEG_AREA_ENRICHMENT      } from './modules/deg_area_enrichment'
include { MASIGPRO_EXPRESSION      } from './modules/masigpro'
include { MASIGPRO_ABUNDANCE       } from './modules/masigpro'
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
params.condition_order    = []      // ordering of conditions; without groups: [test, ref]
params.cutoff_celltypes   = 0
params.delineation_obs_key = null   // null → no delineation extraction
params.areas              = []      // empty → scripts fall back to WHOLE_SAMPLE / ALL_SPOTS
params.masked_celltypes   = ""
params.manual_celltypes   = []
params.sample_metadata    = []      // [{sample, condition, time, group}, ...]
// Time-course design (see deriveContrasts below and masigpro_utils.R::build_design)
params.baseline_group     = null    // group that is the shared time zero (e.g. healthy / SHAM)
params.group_order        = []      // series order; the first is the maSigPro reference (e.g. SCI)
params.contrasts          = []      // [[test, ref], ...]; empty → derived from the groups
params.run_pairwise       = true    // per-contrast CSIDE / DEG_AREA / CT_GENE_EXPR / COMPARE
params.run_timecourse     = true    // maSigPro (needs >= 3 timepoints and a group per sample)
params.gene_list          = null    // legacy single gene list for maSigPro expression
params.masigpro           = [:]     // {areas, gene_lists, degree, Q, alfa, rsq, k, spot_min, ct_min, pseudocount}

// ---------------------------------------------------------------------------
// Derived values from sample_metadata
// ---------------------------------------------------------------------------
def meta            = params.sample_metadata
def samples         = meta.collect { it.sample }
def conditions_map  = meta.collect { "${it.sample}:${it.condition}" }
def timepoints      = meta.collect { it.time }.findAll { it != null }.unique()
def n_timepoints    = timepoints.size()
def withGroup       = meta.count { it.group }
if (withGroup && withGroup != meta.size()) {
    error "sample_metadata: 'group' is set for ${withGroup} of ${meta.size()} samples; set it for all or none"
}
def hasGroups       = withGroup > 0

// Area list (empty → whole-sample fallback inside scripts)
def areaList = params.areas ?: []

// Pairwise contrasts: [[test, ref], ...]
def contrasts = deriveContrasts(meta, params.contrasts ?: [], params.baseline_group,
                                params.group_order ?: [], params.condition_order ?: [])

// maSigPro options
def mopts = params.masigpro ?: [:]
def mo    = { String k, d -> mopts.containsKey(k) && mopts[k] != null ? mopts[k] : d }
def tcAreas   = mo('areas', mo('clusters', areaList))   // 'clusters' = old name
def geneLists = (mo('gene_lists', []) + (params.gene_list ? [params.gene_list] : []))
                    .collect { file(it, checkIfExists: true) }
def tcOpts = [
    baseline_group: params.baseline_group ?: '',
    group_order   : (params.group_order ?: []).join(','),
    degree: mo('degree', 2), Q: mo('Q', 0.05), alfa: mo('alfa', 0.05), rsq: mo('rsq', 0.4), k: mo('k', 9),
    spot_min: mo('spot_min', 0.05), ct_min: mo('ct_min', 0.02), pseudocount: mo('pseudocount', 1e-4)
]
def runTimecourse = params.run_timecourse && n_timepoints >= 3 && hasGroups

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------
workflow {

    log.info "Pairwise contrasts (test_vs_ref): " +
        (params.run_pairwise ? contrasts.collect { "${it[0]}_vs_${it[1]}" }.join(', ') : 'disabled (run_pairwise = false)')
    if (!runTimecourse) {
        log.info "maSigPro skipped: " + (
            !params.run_timecourse ? "run_timecourse = false" :
            n_timepoints < 3       ? "only ${n_timepoints} distinct timepoint(s) ${timepoints} (need >= 3)" :
                                     "sample_metadata has no 'group' per sample")
    }

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
    // 2. Cell2Location: reference per condition, mapping per sample
    //    (memory per job no longer scales with the number of samples)
    // ------------------------------------------------------------------
    def conditionList = params.conditions.tokenize(',').collect { it.trim() }.findAll { it }
    def unknownConds  = meta.collect { it.condition }.unique().findAll { !(it in conditionList) }
    if (unknownConds) {
        error "sample_metadata uses condition(s) ${unknownConds} that are not in params.conditions (${conditionList}); " +
              "those samples would get no cell2location reference"
    }

    CELL2LOC_REFERENCE(
        Channel.fromList(conditionList),
        file(params.single_cell_ref),
        params.ref_label,
        params.masked_celltypes
    )

    // [condition, sample] joined with that condition's reference signatures
    def sampleCondCh = Channel.fromList(meta)
        .map { m -> tuple(m.condition, m.sample) }
    def mapInputCh = sampleCondCh.combine(CELL2LOC_REFERENCE.out.signatures, by: 0)

    CELL2LOC_MAP(
        mapInputCh,
        file(params.spatial_input),
        params.spatial_h5ad_glob,
        params.cutoff_celltypes
    )

    // Gathered once all samples / conditions are done; downstream modules stage
    // these as ./cell2location_map/<sample> and ./data/reference_<cond>
    def cell2locMapCh = CELL2LOC_MAP.out.sample_dir.collect()
    def scRefDataCh   = CELL2LOC_REFERENCE.out.sc_ref_data.collect()

    // ------------------------------------------------------------------
    // 3-6. Pairwise analyses, one task per contrast. Every script below
    //      compares exactly two conditions (condition_order = [test, ref]);
    //      each contrast gets only its own samples and publishes to
    //      <outdir>/contrasts/<test>_vs_<ref>/
    // ------------------------------------------------------------------
    def rctdWeightsCh = Channel.empty()
    if (params.run_pairwise) {
        def contrastCh = Channel.fromList(contrasts.collect { c ->
            def ss = meta.findAll { it.condition in c }
            [ "${c[0]}_vs_${c[1]}".toString(),
              ss.collect { it.sample },
              ss.collect { "${it.sample}:${it.condition}".toString() },
              c ]
        })

        // 3. CSIDE (runs RCTD internally, then CSIDE DE per region)
        CSIDE(
            contrastCh,
            params.ref_level,
            params.ref_label,
            params.masked_celltypes,
            areaList,
            cell2locMapCh,
            EXTRACT_SPATIAL_INPUTS.out.coords,
            delineationCh,
            scRefDataCh
        )
        rctdWeightsCh = CSIDE.out.rctd_weights

        // 4. Per-celltype gene expression & enrichment
        CT_GENE_EXPR_PERCELLTYPE(
            contrastCh.join(CSIDE.out.cside_results),
            params.manual_celltypes,
            cell2locMapCh,
            delineationCh
        )

        CT_GENE_EXPR_ENRICHMENT(
            CT_GENE_EXPR_PERCELLTYPE.out.results
        )

        // 5. DEG by area (independent of the CSIDE chain)
        // NOTE: DEG_areas_pipeline.py reads the raw h5ad (needs full counts + obs),
        //       so it gets spatial_input like CELL2LOC_MAP, not the counts CSVs.
        DEG_AREA_PIPELINE(
            contrastCh,
            file(params.spatial_input),
            params.spatial_h5ad_glob,
            params.ref_level,
            delineationCh,
            areaList
        )

        DEG_AREA_ENRICHMENT(
            DEG_AREA_PIPELINE.out.deg_dir
        )

        // 6. Comparisons (cell2loc vs RCTD/CSIDE)
        COMPARE_CELL2LOC_RCTD(
            contrastCh
                .join(CSIDE.out.rctd_results)
                .join(CSIDE.out.cside_results)
                .join(CT_GENE_EXPR_PERCELLTYPE.out.results)
                .join(CT_GENE_EXPR_PERCELLTYPE.out.cell2loc_map)
                .join(CT_GENE_EXPR_ENRICHMENT.out.results),
            params.ref_level,
            delineationCh,
            areaList
        )
    }

    // ------------------------------------------------------------------
    // 7. Time course (maSigPro): expression per area, and CLR cell-type
    //    abundance from both RCTD and cell2location weights
    // ------------------------------------------------------------------
    if (runTimecourse) {
        def metaCsv = Channel.value(writeMetadataCsv(meta))

        MASIGPRO_EXPRESSION(
            EXTRACT_SPATIAL_INPUTS.out.counts,
            delineationCh,
            metaCsv,
            geneLists,
            tcAreas,
            tcOpts
        )

        // one <sample>_*_weights.csv per sample. RCTD weights are identical
        // across contrasts (per-sample fit), so keep the first copy of each.
        def rctdIn = rctdWeightsCh.flatten()
            .map { f -> [f.name, f] }
            .groupTuple()
            .map { name, fs -> fs[0] }
            .collect()
            .map { fs -> ['rctd', fs] }
        def c2lIn = CELL2LOC_MAP.out.sample_dir
            .map { d -> d.resolve("${d.name}_cell2loc_weights.csv") }
            .collect()
            .map { fs -> ['cell2loc', fs] }
        if (!params.run_pairwise) {
            log.info "maSigPro abundance: RCTD weights come from CSIDE, which is off (run_pairwise = false): cell2location only"
        }

        MASIGPRO_ABUNDANCE(
            c2lIn.mix(rctdIn),
            delineationCh,
            metaCsv,
            tcAreas,
            tcOpts
        )
    }
}

// ---------------------------------------------------------------------------
// Pairwise contrasts [[test, ref], ...]
//   explicit params.contrasts              → used as given
//   groups + baseline_group (time course)  → every condition vs the baseline
//       condition, then at each timepoint every other series vs the reference
//       series (first of group_order), e.g. SCI_10d_vs_healthy, ...,
//       bPAC_10d_vs_SCI_10d, bPAC_30d_vs_SCI_30d
//   otherwise                              → condition_order = [test, ref]
// ---------------------------------------------------------------------------
def deriveContrasts(List meta, List explicit, baseline, List groupOrder, List condOrder) {
    def conds = meta.collect { it.condition }.unique()
    def out = []
    if (explicit) {
        out = explicit.collect { it as List }
    } else if (baseline && meta && meta.every { it.group }) {
        def baseConds = meta.findAll { it.group == baseline }.collect { it.condition }.unique()
        if (baseConds.size() != 1) {
            error "baseline_group '${baseline}' must cover exactly one condition (found ${baseConds}); " +
                  "set params.contrasts explicitly"
        }
        def b = baseConds[0]
        def ordered = condOrder.findAll { it in conds } + conds.findAll { !(it in condOrder) }
        ordered.findAll { it != b }.each { c -> out << [c, b] }

        def series = meta.findAll { it.group != baseline }.collect { it.group }.unique()
        series = groupOrder.findAll { it in series } + series.findAll { !(it in groupOrder) }
        if (series.size() > 1) {
            def ref = series[0]
            meta.findAll { it.group == ref }.collect { it.time }.unique().sort().each { t ->
                def rc = meta.findAll { it.group == ref && it.time == t }.collect { it.condition }.unique()
                series.drop(1).each { g ->
                    def gc = meta.findAll { it.group == g && it.time == t }.collect { it.condition }.unique()
                    if (rc.size() == 1 && gc.size() == 1) out << [gc[0], rc[0]]
                }
            }
        }
    } else {
        if (condOrder.size() != 2) {
            error "Without 'group' + baseline_group in the config, condition_order must be exactly [test, ref] " +
                  "(got ${condOrder}); or list params.contrasts"
        }
        out << condOrder
    }
    out.each { c ->
        if (c.size() != 2 || c[0] == c[1] || !(c[0] in conds) || !(c[1] in conds)) {
            error "Invalid contrast ${c}: need two different conditions from sample_metadata (${conds})"
        }
    }
    return out.unique()
}

// ---------------------------------------------------------------------------
// Helper: write sample_metadata list → CSV file staged for a process
// ---------------------------------------------------------------------------
def writeMetadataCsv(List metadata) {
    def f = file("${workDir}/sample_metadata_${workflow.sessionId}.csv")
    def lines = ["sample,condition,time,group"] +
        metadata.collect { m -> "${m.sample},${m.condition},${m.time},${m.group}" }
    f.text = lines.join("\n") + "\n"
    return f
}
