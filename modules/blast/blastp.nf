nextflow.enable.types = true

// BLASTP of the query sequences against a BLAST database (db: the folder made by MAKEBLASTDB), or against the subject
// sequences, or against themselves when both are null.
// The subject is made into a BLAST database to run on several threads (-num_threads is ignored with -subject).
// Output: tabular format with the query and subject lengths (-outfmt "6 std qlen slen"),
// query.self.blastp (self alignment) or query_vs_subject.blastp (empty when the query is empty).
process BLASTP {
    container 'docker.io/ncbi/blast:2.13.0'
    memory 8.GB
    cpus 4
    time 2.h

    input:
    record(meta: Record, query: Path, subject: Path?, db: Path?, args: String?)

    stage:
    stageAs(db, 'blastdb')

    output:
    record(meta: meta, tsv: file(subject || db ? 'query_vs_subject.blastp' : 'query.self.blastp'))

    script:
    def out = subject || db ? 'query_vs_subject.blastp' : 'query.self.blastp'
    """
    ${db ? '' : "makeblastdb -in '${subject ?: query}' -dbtype prot -out blastdb/db"}
    blastp ${args ?: ''} -num_threads ${task.cpus} -query '${query}' -db blastdb/db -outfmt "6 std qlen slen" -out '${out}'
    """

    stub:
    def out = subject || db ? 'query_vs_subject.blastp' : 'query.self.blastp'
    """
    touch '${out}'
    """
}
