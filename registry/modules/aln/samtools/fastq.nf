nextflow.enable.types = true

process SAMTOOLS_FASTQ {
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    memory 2.GB
    cpus 1
    time 15.m

    input:
    record(meta: Record?, path: Path, args: String?)

    output:
    fastq = record(
        meta: meta,
        reads: file('output.fastq.gz')
    )

    script:
    """
    samtools fastq -@ ${task.cpus} ${args ?: '-n -T "*" -0 -'} ${path} \\
    | bgzip -@ ${task.cpus} \\
    > output.fastq.gz
    """

    stub:
    """
    echo "" | gzip > output.fastq.gz
    """
}
