nextflow.enable.types = true

record MetaPath {
    meta: Record?
    path: Path
}

// BLASTP of the query sequences against the database subject_db (for a self alignment, give the database of the query).
// Output: tabular format with the query and subject lengths (-outfmt "6 std qlen slen"), query_vs_subject.blastp
// (empty when the query is empty).
process BLASTP {
    container 'docker.io/ncbi/blast:2.13.0'
    memory 8.GB
    cpus 4
    time 2.h

    input:
    record(query: MetaPath, subject_db: MetaPath, args: String?)

    stage:
    stageAs(query.path, 'query.fasta')
    stageAs(subject_db.path, 'subject_db')

    output:
    blast6 = record(meta: query.meta, meta_subject: subject_db.meta, path: file('query_vs_subject.blastp'))

    script:
    """
    blastp ${args ?: ''} -num_threads ${task.cpus} -query query.fasta -db subject_db/db -outfmt "6 std qlen slen" -out query_vs_subject.blastp
    """

    stub:
    """
    touch query_vs_subject.blastp
    """
}
