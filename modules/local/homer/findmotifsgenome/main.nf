process HOMER_FINDMOTIFSGENOME {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0fe4a3875b78dce3c66b43fb96489769cc32e55e329e2525d2af09096af2252a/data'
        : 'docker.io/asntech/homer:v5.1'}"

    input:
    tuple val(meta), path(peak_file)
    path fasta
    val genome

    output:
    tuple val(meta), path("${prefix}/"), emit: motifs
    tuple val(meta), path("${prefix}/homerResults.html"), emit: html
    tuple val(meta), path("${prefix}/knownResults.txt"), emit: known, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def size_arg = params.extension_length ? "-size ${params.extension_length}" : "-size 200"
    def genome_ref = fasta ? "${fasta}" : "${genome}"

    """
    findMotifsGenome.pl \\
        ${peak_file} \\
        ${genome_ref} \\
        ${prefix} \\
        ${size_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: \$(echo \$(homer --version 2>&1) | sed 's/^.*v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/homerResults.html
    touch ${prefix}/knownResults.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: 4.11
    END_VERSIONS
    """
}
