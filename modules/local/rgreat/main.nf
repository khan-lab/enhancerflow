process RGREAT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-rgreat:2.4.0--r43hdfd78af_0' :
        'ghcr.io/khan-lab/rgreat:sha-0893d90' }"

    input:
    tuple val(meta), path(super_enhancers)
    val genome

    output:
    tuple val(meta), path("${meta.id}/${meta.id}_great_enrichment.tsv"), emit: go_results
    tuple val(meta), path("${meta.id}/${meta.id}_great_volcano.pdf"), emit: go_volcano, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #Run rGREAT for GO:BP, GO:MF, GO:CC, MSigDB:H 

    run_rgreat.sh \\
        --bed ${super_enhancers} \\
        --genome ${genome} \\
        --outdir ${prefix} \\
        --collection "GO:BP" \\
        --prefix ${prefix}
    
    run_rgreat.sh \\
        --bed ${super_enhancers} \\
        --genome ${genome} \\
        --outdir ${prefix} \\
        --collection "GO:MF" \\
        --prefix ${prefix}_GO_MF

    # Write version information
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    rgreat: \$(Rscript -e 'cat(as.character(utils::packageVersion("rGREAT")))')
    r: \$(Rscript -e 'cat(R.version.string)')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}/${prefix}_great_enrichment.tsv
    touch ${prefix}/${prefix}_great_volcano.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgreat: 2.4.0
    END_VERSIONS
    """
}
