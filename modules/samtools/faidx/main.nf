nextflow.enable.types = true

process SAMTOOLS_FAIDX {
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    memory 2.GB
    cpus 1
    time 15.m

    input:
    record(meta: Record, fasta: Path, args: String?)

    output:
    record(
        meta: meta,
        fai: file("${fasta.name}.fai")
    )

    script:
    """
    samtools faidx ${args ?: ''} '${fasta}'
    """

    stub:
    """
    touch '${fasta.name}.fai'
    """
}
