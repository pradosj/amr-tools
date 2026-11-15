nextflow.enable.types = true

record MetaPath {
    meta: Record?
    path: Path
}


process MINIMAP2_ALIGN_PAF {
    container 'community.wave.seqera.io/library/minimap2_samtools:33bb43c18d22e29c'
    memory 10.GB
    cpus 4
    time 1.h

    input:
    record(meta: Record?, path: Path, ref: MetaPath, args: String?)

    stage:
    stageAs(ref.path, 'ref.fasta')

    output:
    paf = record(meta: meta, meta_ref: ref.meta, path: file('out.paf'))

    script:
    """
    minimap2 ${args ?: ''} -t ${task.cpus} ref.fasta '${path}' > out.paf
    """

    stub:
    """
    touch out.paf
    """
}
