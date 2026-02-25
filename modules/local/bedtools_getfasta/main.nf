process BEDTOOLS_GETFASTA {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::bedtools=2.31.1"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0'
        : 'biocontainers/bedtools:2.31.1--hf5e1c6e_0'}"

    input:
    tuple val(meta), path(bed)
    path fasta

    output:
    tuple val(meta), path("${prefix}.fa"), emit: fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    bedtools getfasta \\
        -fi ${fasta} \\
        -bed ${bed} \\
        ${args} \\
        -fo ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e 's/bedtools v//g')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e 's/bedtools v//g')
    END_VERSIONS
    """
}
