nextflow.enable.types = true

record MetaPath {
    meta: Record?
    path: Path
}


process BWA_MEM {
    container 'community.wave.seqera.io/library/bwa_htslib_samtools:56c9f8d5201889a4'
    errorStrategy task.attempt < 2 ? 'retry' : 'ignore'    // one retry with more resources, then ignore the sample
    memory 8.GB + 8.GB * (task.attempt - 1)
    cpus 4
    time 1.h + 1.h * (task.attempt - 1)

    input:
    record(meta: Record?, read1: Path, read2: Path?, ref: MetaPath, args: String?)    // ref: the folder of BWA_INDEX

    stage:
    stageAs(ref.path, 'index')

    output:
    bam = record(meta: meta, meta_ref: ref.meta, path: file('out.bam'))
    bai = record(meta: meta, meta_ref: ref.meta, path: file('out.bam.bai'))
    outputs = record(meta: meta, meta_ref: ref.meta, bam: file('out.bam'), bai: file('out.bam.bai'))

    script:
    def r2 = read2 ? "${read2}" : ''
    """
    bwa mem ${args ?: ''} -t ${task.cpus} index/index ${read1} ${r2} \\
    | samtools sort -@ ${task.cpus} --write-index -O BAM -o out.bam##idx##out.bam.bai
    """

    stub:
    """
    touch out.bam
    touch out.bam.bai
    """
}
