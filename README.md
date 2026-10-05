# enhancerflow

[](https://github.com/nf-core/enhancerflow#enhancerflow)

> An nf-core-based Nextflow pipeline for identifying and characterising enhancers and super-enhancers from aligned sequencing data.

[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A525.04.0-0DC09D?logo=nextflow&logoColor=white)](https://www.nextflow.io/) [![Docker](https://img.shields.io/badge/Docker-supported-2496ED?logo=docker&logoColor=white)](https://www.docker.com/) [![Conda](https://img.shields.io/badge/Conda-supported-44A833?logo=anaconda&logoColor=white)](https://conda.io/) [![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

The pipeline takes sample-level BAM files and peak calls as input. It calls super-enhancers with ROSE2 and can run downstream visualisation, motif analysis, cis-regulatory circuitry analysis, functional annotation, and comparisons between experimental conditions or time points. A MultiQC report and software-version summary are produced for each run.

## Workflow

The main analysis includes:

- Super-enhancer and typical-enhancer identification from peak files using ROSE2
- Signal visualisation with deepTools
- Motif discovery and enrichment analysis with HOMER, FIMO, and SEA
- Cis-regulatory circuitry analysis with Coltron
- Functional annotation with rGREAT
- Optional cross-condition comparison of enhancer sets and signal tracks
- MultiQC reporting and software-version collection

Individual analysis steps can be disabled with the corresponding `--skip_*` parameters.

## Requirements

- Nextflow 25.04.0 or later
- Docker, Singularity, Apptainer, Conda, or another supported execution environment
- Reference genome files, either through a configured iGenomes genome or explicit paths

Containerised execution is recommended for reproducibility.

## Input samplesheet

The input samplesheet is a comma-separated file. The required columns are `sample`, `peaks`, and `bam`. `control_bam` is optional. `condition` and `timepoint` can be supplied to enable grouped comparisons.

```csv
sample,condition,timepoint,bam,peaks,control_bam
WT_T0_r1,WT,T0,data/WT_T0_r1.bam,data/WT_T0_r1.narrowPeak,data/INPUT_shared.bam
KO_T0_r1,KO,T0,data/KO_T0_r1.bam,data/KO_T0_r1.narrowPeak,data/INPUT_shared.bam
```

Peak files must be in narrowPeak or BED format. BAM files must be coordinate-sorted and accessible from the execution environment.

## Contrastsheet

An optional contrastsheet defines comparisons for differential analysis. It must contain `contrast`, `case`, and `control` columns. Groups are defined by `condition`, or by `condition:timepoint` when time-point information is available.

```csv
contrast,case,control
KO_vs_WT,KO,WT
```

## Running the pipeline

Using an iGenomes reference genome:

```bash
nextflow run nf-core/enhancerflow \
    -profile docker \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38
```

With explicit reference files:

```bash
nextflow run nf-core/enhancerflow \
    -profile docker \
    --input samplesheet.csv \
    --outdir results \
    --fasta /path/to/genome.fa \
    --gtf /path/to/genes.gtf
```

Add an optional contrastsheet with `--contrastsheet contrastsheet.csv`. Use `-resume` to continue a previously interrupted run.

For a quick installation check, run the included test profile:

```bash
nextflow run nf-core/enhancerflow -profile test,docker
```

## Outputs

Results are written to the directory specified by `--outdir`. Depending on the enabled analyses, output includes:

- ROSE2 enhancer and super-enhancer calls and constituent regions
- deepTools signal tracks, matrices, and plots
- Motif, CRC, and functional annotation results
- Cross-condition comparison tables and visualisations
- MultiQC report and pipeline metadata under `pipeline_info/`

## Documentation

See [`docs/usage.md`](docs/usage.md) for additional usage information and [`docs/output.md`](docs/output.md) for output details. Pipeline parameters are defined in [`nextflow_schema.json`](nextflow_schema.json).

## Citations

References for tools used by the pipeline are listed in [`CITATIONS.md`](CITATIONS.md).
