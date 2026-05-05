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
    def flagstat_timeout = params.rose2_flagstat_timeout ?: 60

    """
    # rose2 2.0.1 hard-codes timeout=60 on the samtools flagstat subprocess.run
    # call in rose2/utils.py. Shadow the package on PYTHONPATH with a patched
    # copy so the override works on read-only container filesystems
    # (Singularity/Apptainer) as well as Docker and conda. Children (bamToGFF
    # workers) inherit PYTHONPATH and pick up the same patched module.
    ROSE2_PKG_DIR=\$(python3 -c 'import os, rose2; print(os.path.dirname(rose2.__file__))')
    mkdir -p rose2_override/rose2
    cp -a "\$ROSE2_PKG_DIR"/. rose2_override/rose2/
    sed -i 's/timeout=60)/timeout=${flagstat_timeout})/g' rose2_override/rose2/utils.py
    export PYTHONPATH="\$PWD/rose2_override\${PYTHONPATH:+:\$PYTHONPATH}"

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
