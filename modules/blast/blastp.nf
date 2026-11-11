nextflow.enable.types = true

// Query sequences to align (a FASTA file), and their meta
record BlastSequences {
    meta: Record?
    fasta: Path
}

// BLAST database to align on (the folder made by MAKEBLASTDB with dbtype 'prot': <folder>/db.*), and its meta
record BlastDb {
    meta: Record?
    db: Path
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
    record(query: BlastSequences, subject_db: BlastDb, args: String?)

    stage:
    stageAs(query.fasta, 'query.fasta')
    stageAs(subject_db.db, 'subject_db')

    output:
    record(meta_query: query.meta, meta_subject: subject_db.meta, blast6: file('query_vs_subject.blastp'))

    script:
    """
    blastp ${args ?: ''} -num_threads ${task.cpus} -query query.fasta -db subject_db/db -outfmt "6 std qlen slen" -out query_vs_subject.blastp
    """

    stub:
    """
    touch query_vs_subject.blastp
    """
}
