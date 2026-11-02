nextflow.enable.types = true

process CHOPPER {
    container 'quay.io/biocontainers/chopper:0.12.0--hcdda2d0_0'
    memory 4.GB
    cpus 1
    time 1.h

    input:
    tuple(meta: Map, reads: Path)
    args: String?

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
