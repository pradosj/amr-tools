nextflow.enable.types = true

process FLYE {
    container 'quay.io/biocontainers/flye:2.9.6--py311h2de2dd3_0'
    memory 20.GB
    cpus 8
    time 4.h

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
