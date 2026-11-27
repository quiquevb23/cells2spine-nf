// One task per pairwise contrast, condition_order = [test, ref].
process CT_GENE_EXPR_PERCELLTYPE {
    tag "${cid}"
    label 'process_medium'
    container params.container_rpy2
    publishDir "${params.outdir}/contrasts/${cid}", mode: 'copy', pattern: 'CT_Gene_expr'

    input:
    tuple val(cid), val(samples), val(conditions_map), val(condition_order), path(cside_results)  // "CSIDE" dir of this contrast
    val  manual_celltypes
    path cell2loc_map, stageAs: 'c2l_in/*'   // per-sample dirs from CELL2LOC_MAP
    path delineation_dir      // "delineation" directory from EXTRACT_SPATIAL_INPUTS, or []

    output:
    tuple val(cid), path("CT_Gene_expr"),      emit: results
    // this contrast's cell2location_map/<sample> (upstream files symlinked + the
    // gene_expr_ct_mean/ written here), read by COMPARE_CELL2LOC_RCTD
    tuple val(cid), path("cell2location_map"), emit: cell2loc_map

    script:
    def delin_arg = delineation_dir ? "--delineation_dir ${delineation_dir}" : "--delineation_dir NO_DELINEATION"
    // ct_gene_expr.py writes ./cell2location_map/<sample>/gene_expr_ct_mean/. Its content
    // depends on this contrast's CSIDE gene filter, so build a local cell2location_map of
    // real per-sample dirs (contents symlinked) instead of writing through the staged
    // links into CELL2LOC_MAP's work dirs, which parallel contrasts would share.
    // ct_gene_expr.py only creates CT_Gene_expr/<area> when limma runs (needs >=2 replicates
    // per condition), so make sure the output dir exists even when everything is skipped.
    """
    mkdir -p CT_Gene_expr cell2location_map
    for s in ${samples.join(' ')}; do
        mkdir -p cell2location_map/\$s
        for f in c2l_in/\$s/*; do
            [ "\$(basename \$f)" = gene_expr_ct_mean ] && continue
            ln -s "\$(readlink -f \$f)" cell2location_map/\$s/
        done
    done

    python3 /usr/local/bin/ct_gene_expr.py \\
        --base_dir          . \\
        --out_dir           CT_Gene_expr \\
        --cside_dir         ${cside_results} \\
        --samples           ${samples.join(' ')} \\
        --celltypes         ${manual_celltypes.join(' ')} \\
        --conditions_map    ${conditions_map.join(' ')} \\
        --condition_order   ${condition_order.join(' ')} \\
        ${delin_arg}
    """

    stub:
    """
    mkdir -p CT_Gene_expr/stub
    touch CT_Gene_expr/stub/ct_gene_expr.csv
    for s in ${samples.join(' ')}; do mkdir -p cell2location_map/\$s; done
    """
}
