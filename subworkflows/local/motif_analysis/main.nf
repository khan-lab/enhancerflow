/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: Motif Analysis with HOMER
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HOMER_FINDMOTIFSGENOME      } from '../../../modules/local/homer/findmotifsgenome/main'
include { HOMER_ANNOTATEPEAKS         } from '../../../modules/nf-core/homer/annotatepeaks/main'

workflow MOTIF_ANALYSIS {

    take:
    ch_super_enhancers  // channel: [ meta, super_enhancers_bed ]
    ch_se_constituents  // channel: [ meta, constituent_peaks_bed ]
    ch_fasta            // channel: path(fasta)
    ch_gtf              // channel: path(gtf)
    genome              // val: genome name

    main:
    ch_versions = Channel.empty()

    //
    // Find motifs in SE constituent peaks
    //
    HOMER_FINDMOTIFSGENOME (
        ch_se_constituents,
        ch_fasta,
        genome
    )
    ch_versions = ch_versions.mix(HOMER_FINDMOTIFSGENOME.out.versions.first())

    //
    // Annotate super-enhancers
    //
    HOMER_ANNOTATEPEAKS (
        ch_super_enhancers,
        ch_fasta,
        ch_gtf
    )
    ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS.out.versions.first())

    emit:
    constituents = ch_se_constituents                         // channel: [ meta, bed ]
    motifs       = HOMER_FINDMOTIFSGENOME.out.motifs          // channel: [ meta, motifs ]
    annotations  = HOMER_ANNOTATEPEAKS.out.txt                // channel: [ meta, txt ]
    versions     = ch_versions                                // channel: [ path ]
}
