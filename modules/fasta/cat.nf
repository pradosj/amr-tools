nextflow.enable.types = true

// A FASTA file to concatenate, and the id prefixing the names of its sequences
record IdFasta {
    path: Path
    id: String
}

// Concatenation of FASTA files (in the order of the list), each sequence name prefixed by the id of its file:
// `>name ...` -> `>id:name ...`
// Output: cat.fasta
process FASTA_CAT {
    container 'docker.io/ubuntu:24.04'
    label 'local'
    memory 2.GB
    cpus 1
    time 30.m

    input:
    record(meta: Record, fasta: List<IdFasta>)

    stage:
    stageAs(fasta.collect({ f -> f.path }), 'in/?/*')    // files may share the same name

    output:
    record(meta: meta, fasta: file('cat.fasta'))

    script:
    def cmds = fasta.collect({ f -> "awk -v id='${f.id}' '/^>/ { sub(/^>/, \">\" id \":\") } 1' '${f.path}'" }).join('\n')
    """
    {
    ${cmds}
    } > cat.fasta
    """

    stub:
    """
    touch cat.fasta
    """
}
