process EXTRACT_SPATIAL_INPUTS {
    tag "extraction"
    label 'process_medium'
    container params.container_cell2loc
    publishDir "${params.outdir}/extracted", mode: 'copy'

    input:
    path spatial_input
    val samples
    val delineation_obs_key   // may be null/empty

    output:
    path "counts/*_counts.csv",       emit: counts
    path "coords/*_coords.csv",       emit: coords
    path "delineation/*.csv",         emit: delineation, optional: true

    script:
    def delin_cmd = delineation_obs_key ?
        """
        mkdir -p delineation
        extract_delineation.py \\
            --spatial_input ${spatial_input} \\
            --samples ${samples.join(' ')} \\
            --obs_key ${delineation_obs_key} \\
            --output_dir delineation
        """ : ""
    """
    mkdir -p counts coords
    conversor.py \\
        --spatial_input ${spatial_input} \\
        --output_base_dir . \\
        --samples ${samples.join(' ')}
    ${delin_cmd}
    """
}
