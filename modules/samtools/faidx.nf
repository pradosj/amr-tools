nextflow.enable.types = true

process SAMTOOLS_FAIDX {
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    memory 2.GB
    cpus 1
    time 15.m

    input:
    record(meta: Record?, path: Path, args: String?)

    output:
    fai = record(
        meta: meta,
        path: file("${path.name}.fai")
    )

    script:
    """
    samtools faidx ${args ?: ''} '${path}'
    """

    stub:
    """
    touch '${path.name}.fai'
    """
}
