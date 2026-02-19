process INTERVENE {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::intervene=0.6.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/intervene:0.6.5--pyh5e36f6f_0':
        'biocontainers/intervene:0.6.5--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(beds)
    val overlap_fraction

    output:
    tuple val(meta), path("${prefix}/")           , emit: output_dir
    tuple val(meta), path("${prefix}/Intervene_venn.pdf"), optional: true, emit: venn
    tuple val(meta), path("${prefix}/Intervene_upset.pdf"), optional: true, emit: upset
    tuple val(meta), path("${prefix}/Intervene_pairwise.pdf"), optional: true, emit: heatmap
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def bed_files = beds.collect { it.toString() }.join(' ')

    """
    # Create output directory
    mkdir -p ${prefix}

    # Run intervene for venn diagram (only if ≤6 sets)
    set_count=\$(echo "${bed_files}" | wc -w)
    if [ \$set_count -le 6 ]; then
        intervene venn \\
            -i ${bed_files} \\
            --type genomic \\
            --save-overlaps \\
            -o ${prefix} \\
            ${args}
    fi

    # Run intervene for upset plot
    intervene upset \\
        -i ${bed_files} \\
        --type genomic \\
        --save-overlaps \\
        -o ${prefix} \\
        ${args}

    # Run intervene for pairwise heatmap
    intervene pairwise \\
        -i ${bed_files} \\
        --type genomic \\
        --htype jaccard \\
        -f ${overlap_fraction} \\
        -o ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        intervene: \$(intervene --version | sed 's/intervene //g')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/Intervene_venn.pdf
    touch ${prefix}/Intervene_upset.pdf
    touch ${prefix}/Intervene_pairwise.pdf
    echo '"${task.process}":' > versions.yml
    echo '  intervene: 0.6.5' >> versions.yml
    """
}
