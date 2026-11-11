nextflow.enable.types = true

// The meta of a FASTA file to concatenate: its id, when not null or empty, prefixes the names of its sequences (other fields
// allowed)
record FastaMeta {
	id: String?
}

// A FASTA file to concatenate, and its meta
record IdFasta {
	meta: FastaMeta?    // null: no prefix
	fasta: Path
}

// Concatenation of FASTA files (in the order of the list), each sequence name prefixed by the id of the meta of its
// file when the meta and its id are not null or empty: `>name ...` -> `>id:name ...`. A missing final newline is added.
// Output: cat.fasta
process FASTA_CAT {
    container 'docker.io/ubuntu:24.04'
    label 'local'
    memory 2.GB
    cpus 1
    time 30.m

    input:
    record(meta: Record?, fasta: List<IdFasta>)

    stage:
    stageAs(fasta.collect({ f -> f.fasta }), 'in/?/*')    // files may share the same name

    output:
    record(meta: meta, fasta: file('cat.fasta'))

    script:
    def sep = ':'
    def cmds = fasta.collect({ f -> "awk -v prefix='${f.meta?.id ? "${f.meta.id}${sep}" : ""}' '/^>/ { sub(/^>/, \">\" prefix) } 1' '${f.fasta}'"}).join('\n')
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
