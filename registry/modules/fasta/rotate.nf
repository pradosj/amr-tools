nextflow.enable.types = true

record MetaPath {
    meta: Record?
    path: Path
}

// Rotate every sequence of a genome in a FASTA file according to best hits (highest bit score) found in 
// the given BAM of a library (query) to the FASTA (subject).
process FASTA_ROTATE_BAM {
    container 'registry.gitlab.unige.ch/amr-genomics/rscript:v3'
    label 'rscript'
    memory 8.GB
    cpus 1
    time 30.m

    input:
    record(fasta: MetaPath, bam: MetaPath)

    stage:
    stageAs(fasta.path, 'input.fasta')
    stageAs(bam.path, 'hits.bam')
    stageAs(file("${moduleDir}/assets"), 'assets')

    output:
    fasta = record(meta: fasta.meta, path: file('rotated.fasta'))
    tsv = record(meta: fasta.meta, path: file('rotation.tsv'))

    script:
    """
    #!/usr/bin/env Rscript
    source("assets/lib_io.R")
    source("assets/lib_bam.R")
    source("assets/lib_fasta_rotate.R")
    fa <- fasta_rotate_bam("input.fasta", "hits.bam")
    writeXStringSet(fa, "rotated.fasta")
    write_tsv(as_tibble(attr(fa, "rotation_pos")) |> select(chr = seqnames, pos = start, strand = strand), "rotation.tsv")
    """

    stub:
    """
    touch rotated.fasta
    touch rotation.tsv
    """
}
