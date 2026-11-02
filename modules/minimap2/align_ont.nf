nextflow.enable.types = true

process MINIMAP2_ALIGN_ONT {
    container 'community.wave.seqera.io/library/minimap2_samtools:33bb43c18d22e29c'
    memory 10.GB
    cpus 4
    time 1.h

    input:
    // align the long reads on the assembly
    record(meta: Map, fasta: Path, reads: Path, args: String?)

    output:
    record(
        meta: meta,
        bam: file('out.bam'),
        bai: file('out.bam.bai')
    )

    script:
    """
    minimap2 -x map-ont ${args ?: ''} -t ${task.cpus} ${fasta} ${reads} -a \\
    | samtools sort -@ ${task.cpus} --write-index -O BAM -o out.bam##idx##out.bam.bai
    """

    stub:
    """
    touch out.bam
    touch out.bam.bai
    """
}
