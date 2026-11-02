nextflow.enable.types = true

process MLST {
    container "docker.io/staphb/mlst:2.23.0-2025-02-01"
    memory 4.GB
    cpus 1
    time 30.m

    input:
    record(meta: Map, fasta: Path, args: String?)

    output:
    record(
        meta: meta,
        txt: file('MLST.txt')
    )

    script:
    """
    mlst --threads=${task.cpus} ${args ?: ''} ${fasta} > MLST.txt
    """

    stub:
    """
    touch MLST.txt
    """
}
