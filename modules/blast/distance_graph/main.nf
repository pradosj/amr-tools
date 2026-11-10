nextflow.enable.types = true

// Distance graph of the sequences of a BLAST against themselves (e.g. BLASTN of a FASTA file without subject):
// blast_distance_graph() of R/lib_blast.R, reached through the assets folder of the module.
// Output: blast_distance_graph.rds, an undirected tbl_graph of the sequences, with the distance of each pair (weight:
// errors per aligned bp).
process BLAST_DISTANCE_GRAPH {
    container 'registry.gitlab.unige.ch/amr-genomics/rscript:v3'
    label 'rscript'
    memory 4.GB
    cpus 1
    time 30.m

    input:
    record(meta: Record, blast: Path)

    stage:
    stageAs(blast, 'self.blast')
    stageAs(file("${moduleDir}/assets"), 'assets')    // lib_io.R, lib_blast.R

    output:
    record(
        meta: meta,
        graph: file('blast_distance_graph.rds')
    )

    script:
    """
    #!/usr/bin/env Rscript
    source("assets/lib_io.R")
    source("assets/lib_blast.R")
    g <- blast_distance_graph("self.blast")
    saveRDS(g, "blast_distance_graph.rds")
    """

    stub:
    """
    touch blast_distance_graph.rds
    """
}
