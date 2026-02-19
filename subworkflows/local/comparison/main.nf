/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: Cross-Condition Comparison Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INTERVENE as INTERVENE_INDIVIDUAL   } from '../../../modules/local/intervene/main'
include { INTERVENE as INTERVENE_CONDITIONS   } from '../../../modules/local/intervene/main'
include { DEEPTOOLS_PLOTCORRELATION           } from '../../../modules/nf-core/deeptools/plotcorrelation/main'
include { DEEPTOOLS_COMPUTEMATRIX             } from '../../../modules/nf-core/deeptools/computematrix/main'
include { DEEPTOOLS_PLOTPROFILE               } from '../../../modules/nf-core/deeptools/plotprofile/main'
include { DEEPTOOLS_PLOTHEATMAP               } from '../../../modules/nf-core/deeptools/plotheatmap/main'
include { BEDTOOLS_SORT as SORT_BY_CONDITION  } from '../../../modules/nf-core/bedtools/sort/main'
include { BEDTOOLS_MERGE as MERGE_BY_CONDITION} from '../../../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_SORT as SORT_ALL_SES       } from '../../../modules/nf-core/bedtools/sort/main'
include { BEDTOOLS_MERGE as MERGE_ALL_SES     } from '../../../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_SORT as SORT_ALL_TES       } from '../../../modules/nf-core/bedtools/sort/main'
include { BEDTOOLS_MERGE as MERGE_ALL_TES     } from '../../../modules/nf-core/bedtools/merge/main'

workflow COMPARISON {

    take:
    ch_super_enhancers  // channel: [ meta, bed ] with meta.condition (SEs)
    ch_typical_enhancers // channel: [ meta, bed ] with meta.condition (TEs)
    ch_bigwigs          // channel: [ meta, bigwig ] with meta.condition
    overlap_fraction    // val: minimum overlap fraction

    main:
    ch_versions = Channel.empty()

    //
    // Individual sample comparison with intervene (SEs only, need ≥2 samples)
    //
    ch_individual_ses = ch_super_enhancers
        .map { meta, bed -> bed }
        .collect()
        .filter { beds -> beds.size() >= 2 }
        .map { beds -> [[id: 'all_samples'], beds] }

    INTERVENE_INDIVIDUAL(
        ch_individual_ses,
        overlap_fraction
    )
    ch_versions = ch_versions.mix(INTERVENE_INDIVIDUAL.out.versions.first())

    //
    // Group and concatenate super-enhancers by condition (need ≥2 samples)
    // Gate: only proceed when there are ≥2 SE samples to compare
    //
    ch_se_by_condition = ch_super_enhancers
        .toList()
        .flatMap { items -> items.size() >= 2 ? items : [] }
        .map { meta, bed ->
            def condition = meta.condition ?: 'no_condition'
            [condition, bed]
        }
        .groupTuple()
        .map { condition, beds ->
            def meta = [id: "merged_${condition}", condition: condition]
            [meta, beds]
        }

    // Sort concatenated BEDs
    SORT_BY_CONDITION(
        ch_se_by_condition,
        []  // no genome file needed
    )

    // Merge overlapping regions within each condition
    // Remap meta.id to avoid input/output name collision in bedtools merge
    ch_for_merge_condition = SORT_BY_CONDITION.out.sorted
        .map { meta, bed -> [meta + [id: "${meta.id}_merged"], bed] }

    MERGE_BY_CONDITION(
        ch_for_merge_condition
    )

    //
    // Condition-level comparison with intervene (need ≥2 conditions)
    //
    ch_merged_conditions = MERGE_BY_CONDITION.out.bed
        .map { meta, bed -> bed }
        .collect()
        .filter { beds -> beds.size() >= 2 }
        .map { beds -> [[id: 'condition_comparison'], beds] }

    INTERVENE_CONDITIONS(
        ch_merged_conditions,
        overlap_fraction
    )
    ch_versions = ch_versions.mix(INTERVENE_CONDITIONS.out.versions.first())

    //
    // Collect all bigWigs for multi-sample analysis (need ≥2)
    //
    ch_all_bigwigs = ch_bigwigs
        .map { meta, bigwig -> bigwig }
        .collect()
        .filter { bigwigs -> bigwigs.size() >= 2 }
        .map { bigwigs -> [[id: 'multi_sample'], bigwigs] }

    //
    // Create merged union BED for SEs (need ≥2 for meaningful union)
    //
    ch_all_ses = ch_super_enhancers
        .map { meta, bed -> bed }
        .collect()
        .filter { beds -> beds.size() >= 2 }
        .map { beds -> [[id: 'all_ses_union'], beds] }

    SORT_ALL_SES(
        ch_all_ses,
        []
    )
    ch_for_merge_ses = SORT_ALL_SES.out.sorted
        .map { meta, bed -> [meta + [id: "${meta.id}_merged"], bed] }

    MERGE_ALL_SES(
        ch_for_merge_ses
    )

    //
    // Create merged union BED for TEs (need ≥2 for meaningful union)
    //
    ch_all_tes = ch_typical_enhancers
        .map { meta, bed -> bed }
        .collect()
        .filter { beds -> beds.size() >= 2 }
        .map { beds -> [[id: 'all_tes_union'], beds] }

    SORT_ALL_TES(
        ch_all_tes,
        []
    )
    ch_for_merge_tes = SORT_ALL_TES.out.sorted
        .map { meta, bed -> [meta + [id: "${meta.id}_merged"], bed] }

    MERGE_ALL_TES(
        ch_for_merge_tes
    )

    //
    // deepTools: Compute matrix for all samples at SE + TE regions
    //
    ch_merged_se_bed = MERGE_ALL_SES.out.bed.map { meta, bed -> bed }
    ch_merged_te_bed = MERGE_ALL_TES.out.bed.map { meta, bed -> bed }

    ch_region_beds = ch_merged_se_bed
        .combine(ch_merged_te_bed)
        .map { se_bed, te_bed -> [se_bed, te_bed] }

    ch_for_matrix = ch_all_bigwigs
        .combine(ch_region_beds)

    DEEPTOOLS_COMPUTEMATRIX(
        ch_for_matrix.map { meta, bigwigs, se_bed, te_bed -> [meta, bigwigs] },
        ch_for_matrix.map { meta, bigwigs, se_bed, te_bed -> [se_bed, te_bed] }
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_COMPUTEMATRIX.out.versions.first())

    //
    // deepTools: Correlation from computed matrix
    //
    DEEPTOOLS_PLOTCORRELATION(
        DEEPTOOLS_COMPUTEMATRIX.out.matrix,
        'pearson',  // correlation method
        'heatmap'   // plot type
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_PLOTCORRELATION.out.versions.first())

    //
    // deepTools: Multi-sample heatmap
    //
    DEEPTOOLS_PLOTHEATMAP(
        DEEPTOOLS_COMPUTEMATRIX.out.matrix
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_PLOTHEATMAP.out.versions.first())

    //
    // deepTools: Multi-sample profile plot
    //
    DEEPTOOLS_PLOTPROFILE(
        DEEPTOOLS_COMPUTEMATRIX.out.matrix
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_PLOTPROFILE.out.versions.first())

    emit:
    individual_overlap   = INTERVENE_INDIVIDUAL.out.output_dir      // channel: [ meta, path ]
    individual_venn      = INTERVENE_INDIVIDUAL.out.venn            // channel: [ meta, path ]
    individual_upset     = INTERVENE_INDIVIDUAL.out.upset           // channel: [ meta, path ]
    individual_heatmap   = INTERVENE_INDIVIDUAL.out.heatmap         // channel: [ meta, path ]
    condition_overlap    = INTERVENE_CONDITIONS.out.output_dir      // channel: [ meta, path ]
    condition_venn       = INTERVENE_CONDITIONS.out.venn            // channel: [ meta, path ]
    condition_upset      = INTERVENE_CONDITIONS.out.upset           // channel: [ meta, path ]
    condition_heatmap    = INTERVENE_CONDITIONS.out.heatmap         // channel: [ meta, path ]
    merged_ses           = MERGE_BY_CONDITION.out.bed               // channel: [ meta, bed ]
    correlation_heatmap  = DEEPTOOLS_PLOTCORRELATION.out.pdf        // channel: [ meta, path ]
    multi_heatmap        = DEEPTOOLS_PLOTHEATMAP.out.pdf            // channel: [ meta, path ]
    multi_profile        = DEEPTOOLS_PLOTPROFILE.out.pdf            // channel: [ meta, path ]
    versions             = ch_versions                               // channel: [ path ]
}
