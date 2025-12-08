nextflow.enable.dsl = 2

process VISUALIZE_REPORT {
    publishDir "${params.results_dir}/visualization", mode: 'copy'
    
    input:
    path matrix
    path metadata

    output:
    path "transmission_network.html", emit: network
    path "stats_histogram.png",       emit: hist
    path "stats_heatmap.png",         emit: heatmap
    path "stats_violin.png",          emit: violin

    script:
    """
    python3 $baseDir/scripts/visualize_results.py \\
        --matrix ${matrix} \\
        --metadata ${metadata} \\
        --threshold ${params.snp_threshold}
    """
}

workflow VISUALIZATION {
    take:
    matrix
    metadata

    main:
    VISUALIZE_REPORT(matrix, metadata)
}
