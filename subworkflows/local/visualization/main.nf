/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: Visualization with deepTools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DEEPTOOLS_BAMCOVERAGE   } from '../../../modules/nf-core/deeptools/bamcoverage/main'
include { DEEPTOOLS_COMPUTEMATRIX } from '../../../modules/nf-core/deeptools/computematrix/main'
include { DEEPTOOLS_PLOTPROFILE   } from '../../../modules/nf-core/deeptools/plotprofile/main'
include { DEEPTOOLS_PLOTHEATMAP   } from '../../../modules/nf-core/deeptools/plotheatmap/main'

workflow VISUALIZATION {

    take:
    ch_bams   // channel: [ meta, bam, bai ]
    ch_ses    // channel: [ meta, ses_bed ]
    ch_tes    // channel: [ meta, tes_bed ]

    main:
    ch_versions = Channel.empty()

    //
    // Convert BAMs to bigWig signal tracks
    //
    DEEPTOOLS_BAMCOVERAGE (
        ch_bams,
        [],        // fasta (optional)
        [],        // fasta_fai (optional)
        [[:], []]  // blacklist (optional)
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE.out.versions.first())
    ch_bigwigs = DEEPTOOLS_BAMCOVERAGE.out.bigwig

    //
    // Combine SE and TE BEDs per sample, then join with bigWigs
    //
    ch_region_beds = ch_ses
        .join(ch_tes, by: [0])
        .map { meta, ses_bed, tes_bed -> [meta, [ses_bed, tes_bed]] }

    ch_for_matrix = ch_bigwigs
        .join(ch_region_beds, by: [0])

    //
    // Compute matrix using bigWig files over SE + TE regions
    //
    DEEPTOOLS_COMPUTEMATRIX (
        ch_for_matrix.map { meta, bigwig, beds -> [meta, bigwig] },
        ch_for_matrix.map { meta, bigwig, beds -> beds }
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_COMPUTEMATRIX.out.versions.first())

    //
    // Generate profile plot
    //
    DEEPTOOLS_PLOTPROFILE (
        DEEPTOOLS_COMPUTEMATRIX.out.matrix
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_PLOTPROFILE.out.versions.first())

    //
    // Generate heatmap
    //
    DEEPTOOLS_PLOTHEATMAP (
        DEEPTOOLS_COMPUTEMATRIX.out.matrix
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_PLOTHEATMAP.out.versions.first())

    emit:
    bigwigs      = ch_bigwigs                                   // channel: [ meta, bigwig ]
    matrices     = DEEPTOOLS_COMPUTEMATRIX.out.matrix           // channel: [ meta, mat.gz ]
    profiles     = DEEPTOOLS_PLOTPROFILE.out.pdf                // channel: [ meta, pdf ]
    heatmaps     = DEEPTOOLS_PLOTHEATMAP.out.pdf                // channel: [ meta, pdf ]
    versions     = ch_versions                                  // channel: [ path ]
}
