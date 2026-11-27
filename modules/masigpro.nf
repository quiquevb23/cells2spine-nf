// modules/masigpro.nf
//
// Time-course analyses (paper: maSigPro on area pseudobulk expression and on
// CLR cell-type abundance). Both scripts auto-detect the design from the
// sample metadata CSV (sample, condition, time, group) and exit cleanly with
// SKIPPED.txt when there are < 3 distinct timepoints.

process MASIGPRO_EXPRESSION {
    tag "expression"
    label 'process_medium'
    container params.container_enrichment
    publishDir "${params.outdir}/timecourse", mode: 'copy'

    input:
    path counts_dir              // "counts" dir from EXTRACT_SPATIAL_INPUTS (<sample>_counts.csv)
    path delineation_dir         // "delineation" dir, or [] when not configured
    path sample_metadata         // CSV: sample,condition,time,group
    path gene_lists, stageAs: 'gene_lists/*'   // 0..n gene list files
    val  areas                   // list of delineation labels, or [] for whole sample
    val  opts                    // map: baseline_group, group_order, degree, Q, alfa, rsq, k

    output:
    path "masigpro_expression", emit: results

    script:
    def gl        = [gene_lists].flatten().findAll { it }
    def gene_arg  = gl    ? "--gene_list ${gl.join(',')}"           : "--gene_list NO_FILE"
    def area_arg  = areas ? "--areas ${areas.join(',')}"            : "--areas ALL_SPOTS"
    def delin_arg = delineation_dir ? "--delineation_dir ${delineation_dir}" : "--delineation_dir NO_DELINEATION"
    """
    Rscript /usr/local/bin/masigpro_expression.R \\
        --counts_dir      ${counts_dir} \\
        ${delin_arg} \\
        --sample_metadata ${sample_metadata} \\
        --baseline_group  "${opts.baseline_group}" \\
        --group_order     "${opts.group_order}" \\
        ${area_arg} \\
        ${gene_arg} \\
        --degree ${opts.degree} --Q ${opts.Q} --alfa ${opts.alfa} --rsq ${opts.rsq} --k ${opts.k} \\
        --outdir masigpro_expression
    """

    stub:
    """
    mkdir -p masigpro_expression
    touch masigpro_expression/summary.tsv
    """
}

process MASIGPRO_ABUNDANCE {
    tag "${method}"
    label 'process_low'
    container params.container_enrichment
    publishDir "${params.outdir}/timecourse", mode: 'copy'

    input:
    // method = rctd | cell2loc; weights = flat list of <sample>_RCTD_weights.csv
    // or <sample>_cell2loc_weights.csv (one per sample)
    tuple val(method), path(weights, stageAs: 'weights/*')
    path delineation_dir
    path sample_metadata
    val  areas
    val  opts                    // map: baseline_group, group_order, degree, Q, alfa, rsq, spot_min, ct_min, pseudocount

    output:
    path "masigpro_abundance/${method}", emit: results

    script:
    def pattern   = method == 'rctd' ? '{sample}_RCTD_weights.csv' : '{sample}_cell2loc_weights.csv'
    def area_arg  = areas ? "--areas ${areas.join(',')}" : "--areas ALL_SPOTS"
    def delin_arg = delineation_dir ? "--delineation_dir ${delineation_dir}" : "--delineation_dir NO_DELINEATION"
    """
    Rscript /usr/local/bin/masigpro_abundance.R \\
        --method          ${method} \\
        --weights_dir     weights \\
        --weights_pattern "${pattern}" \\
        ${delin_arg} \\
        --sample_metadata ${sample_metadata} \\
        --baseline_group  "${opts.baseline_group}" \\
        --group_order     "${opts.group_order}" \\
        ${area_arg} \\
        --spot_min ${opts.spot_min} --ct_min ${opts.ct_min} --pseudocount ${opts.pseudocount} \\
        --degree ${opts.degree} --Q ${opts.Q} --alfa ${opts.alfa} --rsq ${opts.rsq} \\
        --outdir masigpro_abundance
    """

    stub:
    """
    mkdir -p masigpro_abundance/${method}
    touch masigpro_abundance/${method}/summary.tsv
    """
}
