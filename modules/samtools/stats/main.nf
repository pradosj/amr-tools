nextflow.enable.types = true

process SAMTOOLS_STATS {
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    memory 12.GB
    cpus 1
    time 30.m

    input:
    record(meta: Map, bam: Path, args: String?)

    output:
    record(
        meta: meta,
        stats: file("${bam.name}.stats")
    )

    script:
    """
    samtools stats -@ ${task.cpus} ${args ?: ''} ${bam} > ${bam.name}.stats
    """

    stub:
    """
    touch ${bam.name}.stats
    """
}
