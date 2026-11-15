nextflow.enable.types = true

process SAMTOOLS_STATS {
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    memory 12.GB
    cpus 1
    time 30.m

    input:
    record(meta: Record?, path: Path, args: String?)

    output:
    stats = record(
        meta: meta,
        path: file("${path.name}.stats")
    )

    script:
    """
    samtools stats -@ ${task.cpus} ${args ?: ''} ${path} > ${path.name}.stats
    """

    stub:
    """
    touch ${path.name}.stats
    """
}
