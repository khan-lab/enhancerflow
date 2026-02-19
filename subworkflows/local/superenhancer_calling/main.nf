/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: Super-Enhancer Calling with ROSE2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ROSE2          } from '../../../modules/local/rose2/main'
include { ROSE_TO_BEDS          } from '../../../modules/local/rose_to_beds/main'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'

workflow SUPERENHANCER_CALLING {

    take:
    ch_samples  // channel: [ meta, peaks, bam, control_bam ]
    //ch_gtf      // channel: path(gtf)
    genome      // val: genome name (e.g., 'hg38', 'mm10')

    main:
    ch_versions = Channel.empty()

//
    // Prepare all BAMs for indexing (both sample and control)
    //
    ch_all_bams = ch_samples
        .flatMap { meta, peaks, bam, control_bam ->
            def bams_to_index = [[meta, bam]]
            if (control_bam && control_bam.name != 'NO_FILE' && !control_bam.toString().isEmpty()) {
                // Create a new meta for control BAM
                def control_meta = meta.clone()
                control_meta.id = "${meta.id}_control"
                bams_to_index.add([control_meta, control_bam])
            }
            return bams_to_index
        }

    //
    // Index all BAM files
    //
    SAMTOOLS_INDEX (
        ch_all_bams
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    //
    // Prepare channel for ROSE2 by joining BAMs with their indices
    //
    ch_bam_bai = SAMTOOLS_INDEX.out.bai
        .map { meta, bai ->
            // Remove _control suffix if present to match original meta
            def original_id = meta.id.replaceAll('_control$', '')
            [original_id, meta.id.endsWith('_control'), bai]
        }

    ch_for_rose2 = ch_samples
        .map { meta, peaks, bam, control_bam ->
            [meta.id, meta, peaks, bam, control_bam]
        }
        .combine(ch_bam_bai, by: 0)
        .map { id, meta, peaks, bam, control_bam, is_control, bai ->
            [id, meta, peaks, bam, control_bam, is_control, bai]
        }
        .groupTuple(by: 0)
        .map { id, meta_list, peaks_list, bam_list, control_bam_list, is_control_list, bai_list ->
            def meta = meta_list[0]
            def peaks = peaks_list[0]
            def bam = bam_list[0]
            def control_bam = control_bam_list[0]

            // Find sample and control BAI files
            def bam_bai = null
            def control_bai = null

            is_control_list.eachWithIndex { is_ctrl, idx ->
                if (is_ctrl) {
                    control_bai = bai_list[idx]
                } else {
                    bam_bai = bai_list[idx]
                }
            }

            // Handle cases where control_bam is empty or null
            if (!control_bam || control_bam.name == 'NO_FILE' || control_bam.toString().isEmpty()) {
                control_bam = file('NO_FILE')
                control_bai = file('NO_FILE')
            }

            [meta, peaks, bam, bam_bai, control_bam, control_bai ?: file('NO_FILE')]
        }

    //
    // Run ROSE2
    //
    ch_all_enhancers   = Channel.empty()
    ch_super_enhancers = Channel.empty()
    ch_plots           = Channel.empty()
    ch_ses              = Channel.empty()
    ch_tes              = Channel.empty()
    ch_ses_constituents = Channel.empty()
    ch_tes_constituents = Channel.empty()

    if (!params.skip_rose2) {
        ROSE2 (
            ch_for_rose2,
            genome
        )
        ch_versions        = ch_versions.mix(ROSE2.out.versions.first())
        ch_all_enhancers   = ROSE2.out.all_enhancers
        ch_super_enhancers = ROSE2.out.super_enhancers
        ch_plots           = ROSE2.out.plot

        //
        // Convert ROSE2 table to SE/TE BED files + constituent peaks
        //
        ch_peaks = ch_samples.map { meta, peaks, bam, control_bam -> [meta, peaks] }
        ch_for_beds = ROSE2.out.all_enhancers.join(ch_peaks, by: [0])

        ROSE_TO_BEDS (
            ch_for_beds
        )
        ch_versions         = ch_versions.mix(ROSE_TO_BEDS.out.versions.first())
        ch_ses              = ROSE_TO_BEDS.out.ses
        ch_tes              = ROSE_TO_BEDS.out.tes
        ch_ses_constituents = ROSE_TO_BEDS.out.ses_constituents
        ch_tes_constituents = ROSE_TO_BEDS.out.tes_constituents
    }

    emit:
    all_enhancers   = ch_all_enhancers   // channel: [ meta, path(txt) ]
    super_enhancers = ch_super_enhancers // channel: [ meta, path(txt) ]
    plots           = ch_plots            // channel: [ meta, path(png) ]
    ses              = ch_ses              // channel: [ meta, path(bed) ]
    tes              = ch_tes              // channel: [ meta, path(bed) ]
    ses_constituents = ch_ses_constituents // channel: [ meta, path(bed) ]
    tes_constituents = ch_tes_constituents // channel: [ meta, path(bed) ]
    bam_bai         = ch_bam_bai          // channel: [ id, is_control, bai ]
    versions        = ch_versions         // channel: [ path ]

}
