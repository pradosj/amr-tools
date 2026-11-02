nextflow.enable.types = true

process BWA_MEM {
    container 'community.wave.seqera.io/library/bwa_htslib_samtools:56c9f8d5201889a4'
    errorStrategy task.attempt < 2 ? 'retry' : 'ignore'    // one retry with more resources, then ignore the sample
    memory 8.GB + 8.GB * (task.attempt - 1)
    cpus 4
    time 1.h + 1.h * (task.attempt - 1)

    input:
    // index: folder produced by BWA_INDEX; reads, read2: short reads
    record(meta: Map, index: Path, reads: Path, read2: Path?, args: String?)

    output:
    record(
        meta: meta,
        bam: file('out.bam'),
        bai: file('out.bam.bai')
    )

    script:
    def r2 = read2 ? "${read2}" : ''
    """
    bwa mem ${args ?: ''} -t ${task.cpus} ${index}/index ${reads} ${r2} \\
    | samtools sort -@ ${task.cpus} --write-index -O BAM -o out.bam##idx##out.bam.bai
    """

    stub:
    """
    touch out.bam
    touch out.bam.bai
    """
}
