#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.fhir_dir        = "$baseDir/data/JSON"
params.reference       = "$baseDir/data/H37Rv.fasta"
params.results_dir     = "$baseDir/results"

log.info """
    TB Phylogeny & Visualization Pipeline (Testing)
"""

include { PHYLO_ANALYSIS } from './workflows/pyhlo.nf'
include { VISUALIZATION }  from './workflows/visualization.nf'

workflow {
    fhir_ch = Channel.fromPath("${params.fhir_dir}/*.json", checkIfExists: true)
    ref_ch = Channel.fromPath(params.reference, checkIfExists: true).first()

    PHYLO_ANALYSIS(fhir_ch, ref_ch)
    VISUALIZATION(PHYLO_ANALYSIS.out.matrix, PHYLO_ANALYSIS.out.metadata)
}