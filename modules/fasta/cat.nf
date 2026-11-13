nextflow.enable.types = true

record IdMeta {
	id: String?
}

record IdPath {
	meta: IdMeta?
	path: Path
}

// Concatenation of FASTA files (in the order of the list), each sequence name prefixed by the id of the meta of its
// file when the meta and its id are not null or empty: `>name ...` -> `>id:name ...`. A missing final newline is added.
// Output: cat.fasta, and cat.tsv, the new id of each sequence (columns meta_id (empty when null), old_seq_id,
// new_seq_id; a sequence id is the first word of its header)
process FASTA_CAT {
    container 'docker.io/ubuntu:24.04'
    label 'local'
    memory 2.GB
    cpus 1
    time 30.m

    input:
    record(meta: Record?, genomes: List<IdPath>)

    stage:
    stageAs(genomes.collect({ f -> f.path }), 'in/?/*')    // files may share the same name

    output:
    fasta = record(meta: meta, path: file('cat.fasta'))
    tsv = record(meta: meta, path: file('cat.tsv'))

    script:
    def sep = ':'
    def cmds = genomes.collect({ f -> """
        awk -v id='${f.meta?.id ?: ''}' \
            -v prefix='${f.meta?.id ? "${f.meta.id}${sep}" : ''}' \
            '/^>/ { old = substr(\$1, 2); sub(/^>/, ">" prefix); print id "\\t" old "\\t" prefix old >> "cat.tsv" } 1' \
            '${f.path}'
    """}).join('\n')
    """
    printf 'meta_id\\told_seq_id\\tnew_seq_id\\n' > cat.tsv
    {
    ${cmds}
    } > cat.fasta
    """

    stub:
    """
    touch cat.fasta cat.tsv
    """
}
