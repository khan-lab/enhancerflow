process ROSE2 {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    
    container 'ghcr.io/khan-lab/rose:2.0.1'

    input:
    tuple val(meta), path(peaks), path(bam), path(bam_index), path(control_bam), path(control_index)
    val genome

    output:
    tuple val(meta), path("*/*_AllStitched.table.txt")  , emit: all_enhancers
    tuple val(meta), path("*/*_SuperStitched.table.txt"), emit: super_enhancers
    tuple val(meta), path("*/*_Plot_points.png")         , emit: plot
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def control = control_bam && control_bam.name != 'NO_FILE' ? "-c ${control_bam}" : ""
    def stitch = params.stitch_distance ?: 12500
    def tss = params.tss_exclusion ?: 2500
    def custom_genome = params.custom_genome ? "--custom ${params.custom_genome}" : ""

    """
    rose2 main -g ${genome.toString().toUpperCase()} \\
        -i ${peaks} \\
        -r ${bam} \\
         ${control} \\
        -o ${prefix} \\
        --tss ${tss} \\
        --stitch ${stitch} \\
        ${custom_genome} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rose2: \$(echo "1.0")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/${prefix}_AllStitched.table.txt
    touch ${prefix}/${prefix}_SuperStitched.table.txt
    touch ${prefix}/${prefix}_Plot_points.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rose2: 1.0
    END_VERSIONS
    """
}
