nextflow.enable.types = true

process BWA_INDEX {
    container 'community.wave.seqera.io/library/bwa_htslib_samtools:56c9f8d5201889a4'
    memory 5.GB
    cpus 3
    time 1.h

    input:
    record(meta: Record?, path: Path, args: String?)

    output:
    bwa_index = record(meta: meta, path: file('index'))

    script:
    """
    mkdir -p "index" && bwa index ${args ?: ''} -p "index/index" ${path}
    """

    stub:
    """
    mkdir -p index
    """
}
