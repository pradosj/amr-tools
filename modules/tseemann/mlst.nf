nextflow.enable.types = true

process MLST {
    container "docker.io/staphb/mlst:2.23.0-2025-02-01"
    memory 4.GB
    cpus 1
    time 30.m

    input:
    record(meta: Record?, path: Path, args: String?)

    output:
    txt = record(meta: meta, path: file('MLST.txt'))

    script:
    """
    mlst --threads=${task.cpus} ${args ?: ''} ${path} > MLST.txt
    """

    stub:
    """
    touch MLST.txt
    """
}
