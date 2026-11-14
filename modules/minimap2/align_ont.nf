nextflow.enable.types = true

record MetaPath {
    meta: Record?
    path: Path
}


process MINIMAP2_ALIGN_ONT {
    container 'community.wave.seqera.io/library/minimap2_samtools:33bb43c18d22e29c'
    memory 10.GB
    cpus 4
    time 1.h

    input:
    record(meta: Record?, path: Path, ref: MetaPath, args: String?)

    stage:
    stageAs(ref.path, 'ref.fasta')

    output:
    bam = record(meta: meta, meta_ref: ref.meta, path: file('out.bam'))
    bai = record(meta: meta, meta_ref: ref.meta, path: file('out.bam.bai'))
    outputs = record(meta: meta, meta_ref: ref.meta, bam: file('out.bam'), bai: file('out.bam.bai'))

    script:
    """
    minimap2 -x map-ont ${args ?: ''} -t ${task.cpus} ref.fasta ${path} -a \\
    | samtools sort -@ ${task.cpus} --write-index -O BAM -o out.bam##idx##out.bam.bai
    """

    stub:
    """
    touch out.bam
    touch out.bam.bai
    """
}
