nextflow.enable.types = true

process BWA_INDEX {
    container 'community.wave.seqera.io/library/bwa_htslib_samtools:56c9f8d5201889a4'
    memory 5.GB
    cpus 3
    time 1.h

    input:
    record(meta: Record, fasta: Path, args: String?)

    output:
    record(
        meta: meta,
        index: file('index')
    )

    script:
    """
    mkdir -p "index" && bwa index ${args ?: ''} -p "index/index" ${fasta}
    """

    stub:
    """
    mkdir -p index
    """
}
