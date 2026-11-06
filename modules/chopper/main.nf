nextflow.enable.types = true

process CHOPPER {
    container 'quay.io/biocontainers/chopper:0.12.0--hcdda2d0_0'
    errorStrategy task.attempt < 2 ? 'retry' : 'ignore'    // one retry with more resources, then ignore the sample
    memory 4.GB + 12.GB * (task.attempt - 1)
    cpus 1
    time 1.h + 1.h * (task.attempt - 1)

    input:
    record(meta: Record, reads: Path, args: String?)

    output:
    record(
        meta: meta,
        reads: file('chopped_reads.fastq.gz')
    )

    script:
    """
    chopper ${args ?: '--quality 15 --minlength 200'} \\
        --input ${reads} \\
        --threads ${task.cpus} \\
    | gzip --fast \\
    > chopped_reads.fastq.gz
    """

    stub:
    """
    echo "" | gzip > chopped_reads.fastq.gz
    """
}
