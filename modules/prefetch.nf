process PREFETCH {
    tag "$srr"

    input:
    tuple val(srr), val(sample_id), val(type)

    output:
    tuple val(srr), val(sample_id), val(type), path("${srr}/*"), emit: sra_dir

    script:
    """
    prefetch ${srr}
    """
}
