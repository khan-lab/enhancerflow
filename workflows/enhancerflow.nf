/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap            } from 'plugin/nf-schema'
include { paramsSummaryMultiqc        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText      } from '../subworkflows/local/utils_nfcore_enhancerflow_pipeline'
include { SAMTOOLS_INDEX              } from '../modules/nf-core/samtools/index/main'

include { SUPERENHANCER_CALLING       } from '../subworkflows/local/superenhancer_calling/main'
include { VISUALIZATION               } from '../subworkflows/local/visualization/main'
include { MOTIF_ANALYSIS              } from '../subworkflows/local/motif_analysis/main'
include { CRC_ANALYSIS                } from '../subworkflows/local/crc_analysis/main'
include { FUNCTIONAL_ANNOTATION       } from '../subworkflows/local/functional_annotation/main'
include { COMPARISON                  } from '../subworkflows/local/comparison/main'
include { GENERATE_DESIGN_SUMMARY     } from '../modules/local/generate_design_summary/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ENHANCERFLOW {

    take:
    ch_samplesheet   // channel: [ meta, peaks, bam, control_bam ]
    ch_contrastsheet // channel: [ meta, case_ids, control_ids ] or empty

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Generate design summary if contrastsheet provided
    //
    if (params.contrastsheet) {
        GENERATE_DESIGN_SUMMARY (
            ch_contrastsheet
        )
        ch_versions = ch_versions.mix(GENERATE_DESIGN_SUMMARY.out.versions)
    }

    //
    // Prepare genome files
    //
    ch_fasta = params.fasta ? Channel.fromPath(params.fasta) : Channel.value([])
    ch_gtf   = params.gtf   ? Channel.fromPath(params.gtf)   : Channel.value([])

    //
    // SUBWORKFLOW: Super-Enhancer Calling with ROSE2
    //
    SUPERENHANCER_CALLING (
        ch_samplesheet,
        //ch_gtf,
        params.genome
    )
    ch_versions = ch_versions.mix(SUPERENHANCER_CALLING.out.versions)

    //
    // SUBWORKFLOW: Visualization with deepTools
    //
    // SUBWORKFLOW: Visualization with deepTools
    //
    if (!params.skip_visualization) {
        // Index BAMs for visualization
        ch_bams_for_viz = ch_samplesheet
            .map { meta, peaks, bam, control_bam -> [meta, bam] }

        SAMTOOLS_INDEX (
            ch_bams_for_viz
        )

        // Combine BAMs with their indices
        ch_bams_indexed = ch_bams_for_viz
            .join(SAMTOOLS_INDEX.out.bai, by: 0)
            .map { meta, bam, bai -> [meta, bam, bai] }

        VISUALIZATION (
            ch_bams_indexed,
            SUPERENHANCER_CALLING.out.ses,
            SUPERENHANCER_CALLING.out.tes
        )
        ch_versions = ch_versions.mix(VISUALIZATION.out.versions)
    }

    //
    // SUBWORKFLOW: Motif Analysis with HOMER
    //
    if (!params.skip_motifs) {
        MOTIF_ANALYSIS (
            SUPERENHANCER_CALLING.out.ses,
            SUPERENHANCER_CALLING.out.ses_constituents,
            ch_fasta,
            ch_gtf,
            params.genome
        )
        ch_versions = ch_versions.mix(MOTIF_ANALYSIS.out.versions)
    }

    //
    // SUBWORKFLOW: CRC Analysis with Coltron
    //
    if (!params.skip_crc && !params.skip_motifs && params.fasta) {
        CRC_ANALYSIS (
            SUPERENHANCER_CALLING.out.all_enhancers,
            MOTIF_ANALYSIS.out.constituents,
            ch_fasta,
            params.genome
        )
        ch_versions = ch_versions.mix(CRC_ANALYSIS.out.versions)
    }

    //
    // SUBWORKFLOW: Functional Annotation with rGREAT
    //
    if (!params.skip_great) {
        FUNCTIONAL_ANNOTATION (
            SUPERENHANCER_CALLING.out.ses,
            params.genome
        )
        ch_versions = ch_versions.mix(FUNCTIONAL_ANNOTATION.out.versions)
    }

    //
    // SUBWORKFLOW: Cross-Condition Comparison
    //
    // Requires BAM files with condition metadata
    if (!params.skip_comparison && !params.skip_visualization) {
        // Filter for samples with condition metadata
        ch_se_with_conditions = SUPERENHANCER_CALLING.out.ses
            .filter { meta, bed -> meta.condition != null }

        ch_te_with_conditions = SUPERENHANCER_CALLING.out.tes
            .filter { meta, bed -> meta.condition != null }

        ch_bw_with_conditions = VISUALIZATION.out.bigwigs
            .filter { meta, bigwig -> meta.condition != null }

        // Count samples with conditions for logging
        ch_se_with_conditions
            .count()
            .subscribe { count ->
                if (count > 1) {
                    log.info "Running cross-condition comparison analysis on ${count} samples with condition metadata..."
                } else {
                    log.info "No samples with condition metadata found. Skipping comparison analysis."
                }
            }

        COMPARISON (
            ch_se_with_conditions,
            ch_te_with_conditions,
            ch_bw_with_conditions,
            params.overlap_fraction
        )
        ch_versions = ch_versions.mix(COMPARISON.out.versions)
    }

    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'enhancerflow_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
