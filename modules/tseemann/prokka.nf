nextflow.enable.types = true

process PROKKA {
    container "docker.io/staphb/prokka:1.14.6"
    memory 12.GB
    cpus 4
    time 1.h

    input:
    record(meta: Map, fasta: Path, args: String?)

    output:
    record(
        meta: meta,
        dir: file('prokka')
    )

    script:
    """
    prokka ${args ?: ''} --outdir 'prokka/' --cpus ${task.cpus} --prefix prokka ${fasta}
    """

    stub:
    """
    mkdir -p prokka
    """
}
