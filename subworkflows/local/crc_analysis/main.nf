/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: Core Regulatory Circuitry Analysis with Coltron
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { COLTRON } from '../../../modules/local/coltron/main'

workflow CRC_ANALYSIS {

    take:
    ch_enhancers     // channel: [ meta, all_enhancers ]
    ch_constituents  // channel: [ meta, constituent_bed ]
    ch_fasta         // channel: path(fasta)
    genome           // val: genome name

    main:
    ch_versions = Channel.empty()

    //
    // Combine enhancer table with constituent peaks
    //
    ch_for_coltron = ch_enhancers.join(ch_constituents, by: [0])

    //
    // Run Coltron CRC analysis
    //
    COLTRON (
        ch_for_coltron,
        ch_fasta,
        genome
    )
    ch_versions = ch_versions.mix(COLTRON.out.versions.first())

    emit:
    output_dir    = COLTRON.out.output_dir     // channel: [ meta, dir ]
    clique_scores = COLTRON.out.clique_scores  // channel: [ meta, txt ]
    degree_table  = COLTRON.out.degree_table   // channel: [ meta, txt ]
    edge_list     = COLTRON.out.edge_list      // channel: [ meta, txt ]
    versions      = ch_versions                // channel: [ path ]
}
