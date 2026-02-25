/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: Motif Analysis with HOMER, FIMO, and SEA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HOMER_FINDMOTIFSGENOME as HOMER_FINDMOTIFSGENOME_SE } from '../../../modules/local/homer/findmotifsgenome/main'
include { HOMER_FINDMOTIFSGENOME as HOMER_FINDMOTIFSGENOME_TE } from '../../../modules/local/homer/findmotifsgenome/main'
include { HOMER_ANNOTATEPEAKS                   } from '../../../modules/nf-core/homer/annotatepeaks/main'
include { BEDTOOLS_GETFASTA as GETFASTA_SE      } from '../../../modules/local/bedtools_getfasta/main'
include { BEDTOOLS_GETFASTA as GETFASTA_TE      } from '../../../modules/local/bedtools_getfasta/main'
include { FIMO as FIMO_SE                       } from '../../../modules/local/fimo/main'
include { FIMO as FIMO_TE                       } from '../../../modules/local/fimo/main'
include { SEA as SEA_SE                         } from '../../../modules/local/sea/main'
include { SEA as SEA_TE                         } from '../../../modules/local/sea/main'

workflow MOTIF_ANALYSIS {

    take:
    ch_super_enhancers  // channel: [ meta, super_enhancers_bed ]
    ch_se_constituents  // channel: [ meta, constituent_peaks_bed ]
    ch_te_constituents  // channel: [ meta, constituent_peaks_bed ]
    ch_fasta            // channel: path(fasta)
    ch_gtf              // channel: path(gtf)
    ch_motif_db         // channel: path(motif_db)
    genome              // val: genome name

    main:
    ch_versions = Channel.empty()

    //
    // Find motifs in SE constituent peaks with HOMER
    //
    HOMER_FINDMOTIFSGENOME_SE (
        ch_se_constituents,
        ch_fasta,
        genome
    )
    ch_versions = ch_versions.mix(HOMER_FINDMOTIFSGENOME_SE.out.versions.first())

    //
    // Find motifs in TE constituent peaks with HOMER
    //
    HOMER_FINDMOTIFSGENOME_TE (
        ch_te_constituents,
        ch_fasta,
        genome
    )
    ch_versions = ch_versions.mix(HOMER_FINDMOTIFSGENOME_TE.out.versions.first())

    //
    // Annotate super-enhancers with HOMER
    //
    HOMER_ANNOTATEPEAKS (
        ch_super_enhancers,
        ch_fasta,
        ch_gtf
    )
    ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS.out.versions.first())

    //
    // Convert SE and TE constituent BEDs to FASTA
    //
    GETFASTA_SE (
        ch_se_constituents,
        ch_fasta
    )
    ch_versions = ch_versions.mix(GETFASTA_SE.out.versions.first())

    GETFASTA_TE (
        ch_te_constituents,
        ch_fasta
    )
    ch_versions = ch_versions.mix(GETFASTA_TE.out.versions.first())

    //
    // FIMO: Scan for individual motif occurrences
    //
    FIMO_SE (
        GETFASTA_SE.out.fasta,
        ch_motif_db
    )
    ch_versions = ch_versions.mix(FIMO_SE.out.versions.first())

    FIMO_TE (
        GETFASTA_TE.out.fasta,
        ch_motif_db
    )
    ch_versions = ch_versions.mix(FIMO_TE.out.versions.first())

    //
    // SEA: Simple Enrichment Analysis
    //
    SEA_SE (
        GETFASTA_SE.out.fasta,
        ch_motif_db
    )
    ch_versions = ch_versions.mix(SEA_SE.out.versions.first())

    SEA_TE (
        GETFASTA_TE.out.fasta,
        ch_motif_db
    )
    ch_versions = ch_versions.mix(SEA_TE.out.versions.first())

    emit:
    constituents = ch_se_constituents                         // channel: [ meta, bed ]
    homer_motifs_se    = HOMER_FINDMOTIFSGENOME_SE.out.motifs        // channel: [ meta, motifs ]
    homer_motifs_te    = HOMER_FINDMOTIFSGENOME_TE.out.motifs        // channel: [ meta, motifs ]
    homer_annotations  = HOMER_ANNOTATEPEAKS.out.txt                // channel: [ meta, txt ]
    se_fimo      = FIMO_SE.out.tsv                            // channel: [ meta, tsv ]
    te_fimo      = FIMO_TE.out.tsv                            // channel: [ meta, tsv ]
    se_sea       = SEA_SE.out.tsv                             // channel: [ meta, tsv ]
    te_sea       = SEA_TE.out.tsv                             // channel: [ meta, tsv ]
    versions     = ch_versions                                // channel: [ path ]
}
