nextflow.enable.types = true

// BLASTP of the query sequences against the subject sequences, or against themselves when subject is null.
// The subject is made into a BLAST database to run on several threads (-num_threads is ignored with -subject).
// Output: tabular format with the query and subject lengths (-outfmt "6 std qlen slen"),
// query.self.blastp (self alignment) or query_vs_subject.blastp (empty when the query is empty).
process BLASTP {
    container 'docker.io/ncbi/blast:2.13.0'
    memory 8.GB
    cpus 4
    time 2.h

    input:
    record(meta: Map, query: Path, subject: Path?, args: String?)

    output:
    record(meta: meta, tsv: file(subject ? 'query_vs_subject.blastp' : 'query.self.blastp'))

    script:
    def out = subject ? 'query_vs_subject.blastp' : 'query.self.blastp'
    """
    makeblastdb -in '${subject ?: query}' -dbtype prot -out blastdb/subject
    blastp ${args ?: ''} -num_threads ${task.cpus} -query '${query}' -db blastdb/subject -outfmt "6 std qlen slen" -out '${out}'
    """

    stub:
    def out = subject ? 'query_vs_subject.blastp' : 'query.self.blastp'
    """
    touch '${out}'
    """
}
