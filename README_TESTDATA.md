# EnhancerFlow minimal test bundle (synthetic)

This bundle is intended for **schema/input parsing tests** and contrasts logic.
It includes synthetic peak files and **placeholder empty BAM/BAI files**.

If your pipeline runs tools that require valid BAMs (samtools/deeptools/etc.),
use these files only in a test profile that **skips BAM-consuming steps**, or
replace BAMs with real small BAMs.

## Files
- `samplesheet.csv`
- `contrastsheet_condition.csv`
- `contrastsheet_longitudinal.csv`
- `data/*.narrowPeak`
- `data/*.bam` and `data/*.bam.bai` (empty placeholders)
- `data/test_genome.fa` (optional reference)
