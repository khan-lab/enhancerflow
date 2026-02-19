process COLTRON {
    tag "${meta.id}"
    label 'process_high'

    container "shin111/coltron:version1.2"

    input:
    tuple val(meta), path(enhancer_table), path(subpeak_bed)
    path fasta
    val genome

    output:
    tuple val(meta), path("${prefix}/"), emit: output_dir
    tuple val(meta), path("${prefix}/*_CLIQUE_SCORES*.txt"), emit: clique_scores, optional: true
    tuple val(meta), path("${prefix}/*_DEGREE_TABLE.txt"), emit: degree_table, optional: true
    tuple val(meta), path("${prefix}/*_EDGE_LIST.txt"), emit: edge_list, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // Get directory containing chromosome FASTA files
    def fasta_dir = fasta.getParent()

    """
    mkdir -p ${prefix}

    coltron \\
        -e ${enhancer_table} \\
        -g ${genome} \\
        -s ${subpeak_bed} \\
        -c ${fasta_dir} \\
        -o ${prefix} \\
        -n ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coltron: \$(echo "1.2")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/${prefix}_CLIQUE_SCORES.txt
    touch ${prefix}/${prefix}_DEGREE_TABLE.txt
    touch ${prefix}/${prefix}_EDGE_LIST.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coltron: 1.2
    END_VERSIONS
    """
}
