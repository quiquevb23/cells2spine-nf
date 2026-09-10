process EXTRACT_SPATIAL_INPUTS {
    tag "extraction"
    label 'process_medium'
    container params.container_cell2loc
    publishDir "${params.outdir}/extracted", mode: 'copy'

    input:
    path spatial_input
    val  samples
    val  delineation_obs_key   // may be null/empty
    val  h5ad_glob

    output:
    path "counts",      emit: counts
    path "coords",      emit: coords
    path "delineation", emit: delineation, optional: true

    script:
    def delin_cmd = delineation_obs_key ?
        """
        mkdir -p delineation
        extract_delineation.py \\
            --spatial_input ${spatial_input} \\
            --samples ${samples.join(' ')} \\
            --obs_key ${delineation_obs_key} \\
            --output_dir delineation \\
            --spatial_h5ad_glob "${h5ad_glob}"
        """ : ""
    """
    mkdir -p counts coords

    conversor.py \\
        --spatial_input ${spatial_input} \\
        --output_base_dir . \\
        --samples ${samples.join(' ')} \\
        --spatial_h5ad_glob "${h5ad_glob}"

    ${delin_cmd}
    """

    stub:
    """
    mkdir -p counts coords delineation
    touch counts/placeholder coords/placeholder delineation/placeholder
    """
}
