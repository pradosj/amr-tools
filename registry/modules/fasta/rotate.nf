nextflow.enable.types = true

record MetaPath {
    meta: Record?
    path: Path
}

// Rotate every sequence of a genome in a FASTA file according to best hits (highest bit score) found in 
// the given BLAST to a reference library.
process FASTA_ROTATE {
    container 'registry.gitlab.unige.ch/amr-genomics/rscript:v3'
    label 'rscript'
    memory 8.GB
    cpus 1
    time 30.m

    input:
    record(fasta: MetaPath, blast6: MetaPath)

    stage:
    stageAs(fasta.path, 'input.fasta')
    stageAs(blast6.path, 'hits.blast6')
    stageAs(file("${moduleDir}/assets"), 'assets')

    output:
    record(meta: fasta.meta, path: file('rotated.fasta'))

    script:
    """
    #!/usr/bin/env Rscript
    source("assets/lib_io.R")
    source("assets/lib_fasta_rotate.R")
    rotate_fasta("input.fasta", "hits.blast6", "rotated.fasta")
    """

    stub:
    """
    touch rotated.fasta
    """
}
