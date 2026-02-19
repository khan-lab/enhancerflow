process EXTRACT_CONSTITUENTS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(all_enhancers), path(peaks)

    output:
    tuple val(meta), path("*_constituent_peaks.bed"), emit: constituents
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3

    import sys

    # Read super-enhancer regions
    se_regions = []
    with open("${all_enhancers}") as f:
        for line in f:
            if line.startswith('REGION_ID') or line.startswith('#'):
                continue
            parts = line.strip().split('\\t')
            if len(parts) >= 3:
                try:
                    chrom, start, end = parts[1].split(':')[0], int(parts[1].split(':')[1].split('-')[0]), int(parts[1].split(':')[1].split('-')[1])
                    se_regions.append((chrom, start, end))
                except:
                    continue

    # Extract peaks overlapping SEs
    with open("${peaks}") as fin, open("${prefix}_constituent_peaks.bed", "w") as fout:
        for line in fin:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\\t')
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])

            # Check overlap with any SE
            for se_chr, se_start, se_end in se_regions:
                if chrom == se_chr and not (end < se_start or start > se_end):
                    fout.write(line)
                    break

    # Versions
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: {sys.version.split()[0]}\\n')
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_constituent_peaks.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9.0
    END_VERSIONS
    """
}
