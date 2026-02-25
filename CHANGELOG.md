# nf-core/enhancerflow: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0dev - [TBD]

### Major Changes

- Pipeline renamed from nf-core/seia to nf-core/enhancerflow
- Added contrastsheet functionality for differential analysis design
- Added FIMO and SEA motif analysis using MEME suite
- Added ROSE_TO_BEDS module for generating SE/TE BED files and constituent peaks

### `Added`

- Optional `timepoint` column in samplesheet for temporal grouping
- Contrastsheet parameter for defining differential analysis contrasts
- Design summary generation for contrasts
- Contrast-specific comparison analysis in COMPARISON subworkflow
- ROSE_TO_BEDS module to convert ROSE2 output to SE/TE BED files with constituent peaks
- BEDTOOLS_GETFASTA module to extract FASTA sequences from BED regions
- FIMO module for scanning sequences for individual motif occurrences
- SEA module for simple enrichment analysis of known motifs
- `--motif_db` parameter for custom motif database (defaults to JASPAR in MEME container)
- DEEPTOOLS_BAMCOVERAGE step to convert BAMs to bigWigs before computeMatrix
- Both SE and TE BED regions passed to computeMatrix

### `Fixed`

- Fixed genome parameter overriding CLI values (use Elvis operator in main.nf)
- Fixed Channel.empty() preventing HOMER and COLTRON from running (changed to Channel.value([]))
- Fixed Channel.fromPath() creating queue channels consumed after first use (changed to Channel.value(file()))
- Fixed HOMER_FINDMOTIFSGENOME/ANNOTATEPEAKS to handle empty fasta/gtf gracefully
- Fixed COMPARISON subworkflow name collision in bedtools merge (sort/merge meta.id remapping)
- Added minimum sample guards in COMPARISON subworkflow (skip with <2 samples)
- Guarded CRC_ANALYSIS to only run when fasta is provided

### `Dependencies`

- MEME suite (FIMO, SEA) via `ghcr.io/khan-lab/meme`
- HOMER v5.1 via `docker.io/asntech/homer:v5.1`

### `Deprecated`

## v1.0.0dev - [date]

Initial release of nf-core/enhancerflow, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
