nextflow.enable.types = true

// Concatenate FASTA files, in the given order (adding the final newline missing in some files)
process FASTA_CONCAT {
    container 'docker.io/ubuntu:24.04'
    label 'local'
    memory 2.GB
    cpus 1
    time 30.m

    input:
    record(meta: Map, fasta: List<Path>)

    stage:
    stageAs(fasta, 'in/?/*')    // files may share the same name

    output:
    record(meta: meta, fasta: file('concat.fasta'))

    script:
    """
    awk 1 ${fasta.collect({ f -> "'${f}'" }).join(' ')} > concat.fasta
    """

    stub:
    """
    touch concat.fasta
    """
}
