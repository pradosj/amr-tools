nextflow.enable.types = true

record MetaPath {
    meta: Record?
    path: Path
}

// Rotation of the sequences of a FASTA file: each sequence starts at its best hit (highest bit score) of the BLAST hits,
// on the strand of the hit (sequences without hit are kept as they are). blast6: hits of reference sequences (query) on
// the sequences to rotate (subject), -outfmt "6 std qlen slen" (see BLASTN).
// Output: rotated.fasta, with the meta of the FASTA file. The R functions are in assets/lib_fasta_rotate.R.
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
