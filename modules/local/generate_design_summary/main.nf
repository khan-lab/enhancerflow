process GENERATE_DESIGN_SUMMARY {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/64/6494f63c04a69da6baf9a70769f6cbb95c5d1b8a7ae99d4c9a3c04f96bbdfbc1/data' :
        'community.wave.seqera.io/library/bash:5.2.37--06a61d9eb4b81d26' }"

    input:
    tuple val(meta), val(case_sample_ids), val(control_sample_ids)

    output:
    path "${prefix}_design_summary.tsv", emit: summary
    path "versions.yml"                 , emit: versions

    script:
    def prefix = meta.id
    def case_samples = case_sample_ids.join(',')
    def control_samples = control_sample_ids.join(',')
    """
    echo -e "contrast\\tcase_group\\tcontrol_group\\tsample_ids_in_case\\tsample_ids_in_control" > ${prefix}_design_summary.tsv
    echo -e "${meta.contrast_id}\\t${meta.case_group}\\t${meta.control_group}\\t${case_samples}\\t${control_samples}" >> ${prefix}_design_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | cut -d' ' -f4)
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch ${prefix}_design_summary.tsv
    echo '"${task.process}":' > versions.yml
    echo '  bash: 5.0.0' >> versions.yml
    """
}
