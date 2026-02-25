process FIMO {
    tag "${meta.id}"
    label 'process_medium'

    container "ghcr.io/khan-lab/meme:sha-3aeef47"

    input:
    tuple val(meta), path(sequences)
    path motif_db

    output:
    tuple val(meta), path("${prefix}_fimo/"),         emit: output_dir
    tuple val(meta), path("${prefix}_fimo/fimo.tsv"), emit: tsv
    tuple val(meta), path("${prefix}_fimo/fimo.html"),emit: html, optional: true
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def db = motif_db ? "${motif_db}" : "/opt/meme/db/motif_databases/JASPAR/JASPAR2024_CORE_vertebrates_non-redundant_v2.meme"
    """
    fimo \\
        --oc ${prefix}_fimo \\
        ${args} \\
        ${db} \\
        ${sequences}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meme: \$(fimo --version 2>&1 | head -1 | sed 's/.* //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_fimo
    touch ${prefix}_fimo/fimo.tsv
    touch ${prefix}_fimo/fimo.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meme: 5.5.5
    END_VERSIONS
    """
}
