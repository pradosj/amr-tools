nextflow.enable.types = true

// Align an assembly (query) to a reference assembly
process MINIMAP2_ALIGN_ASM5 {
    container 'community.wave.seqera.io/library/minimap2_samtools:33bb43c18d22e29c'
    memory 10.GB
    cpus 4
    time 1.h

    input:
    record(meta: Record, ref: Path, query: Path, args: String?)

    stage:
    stageAs(ref, 'ref.fasta')
    stageAs(query, 'asm.fasta')

    output:
    record(meta: meta, bam: file('out.bam'), bai: file('out.bam.bai'))

    script:
    """
    # the same file given as ref and query is staged only once
    [ -e asm.fasta ] || ln -s ref.fasta asm.fasta
    minimap2 -x asm5 ${args ?: ''} -t ${task.cpus} ref.fasta asm.fasta -a \\
        | samtools sort -@ ${task.cpus} --write-index -O BAM -o out.bam##idx##out.bam.bai
    """

    stub:
    """
    [ -e asm.fasta ] || ln -s ref.fasta asm.fasta
    touch out.bam out.bam.bai
    """
}
