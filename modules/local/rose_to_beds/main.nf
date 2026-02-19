process ROSE_TO_BEDS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(all_enhancers), path(peaks)

    output:
    tuple val(meta), path("${prefix}_SEs.bed"),                emit: ses
    tuple val(meta), path("${prefix}_TEs.bed"),                emit: tes
    tuple val(meta), path("${prefix}_SEs_constituents.bed"),   emit: ses_constituents
    tuple val(meta), path("${prefix}_TEs_constituents.bed"),   emit: tes_constituents
    path "versions.yml",                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python ${moduleDir}/bin/rose_to_beds.py \\
        --rose ${all_enhancers} \\
        --peaks ${peaks} \\
        --outdir . \\
        --prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_SEs.bed
    touch ${prefix}_TEs.bed
    touch ${prefix}_SEs_constituents.bed
    touch ${prefix}_TEs_constituents.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9.0
    END_VERSIONS
    """
}
