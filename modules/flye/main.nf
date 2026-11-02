nextflow.enable.types = true

process FLYE {
    container 'quay.io/biocontainers/flye:2.9.6--py311h2de2dd3_0'
    errorStrategy task.attempt < 2 ? 'retry' : 'ignore'    // one retry with more resources, then ignore the sample
    memory 20.GB + 32.GB * (task.attempt - 1)
    cpus 8
    time 4.h + 8.h * (task.attempt - 1)

    input:
    record(meta: Map, reads: Path, args: String?)

    output:
    record(
        meta: meta,
        dir: file('flye'),
        fasta: file('flye/assembly.fasta')
    )

    script:
    """
    flye \\
        ${args ?: ''} \\
        --threads ${task.cpus} \\
        --out-dir ./flye \\
        --nano-hq ${reads}
    """

    stub:
    """
    mkdir -p flye/
    touch flye/assembly.fasta
    """
}
