nextflow.enable.types = true

// BLASTN of the query sequences against the subject sequences, or against themselves when subject is null.
// The subject is made into a BLAST database to run on several threads (-num_threads is ignored with -subject).
// Output: tabular format with the query and subject lengths (-outfmt "6 std qlen slen"),
// query.self.blastn (self alignment) or query_vs_subject.blastn (empty when the query is empty).
process BLASTN {
    container 'docker.io/ncbi/blast:2.13.0'
    memory 8.GB
    cpus 4
    time 2.h

    input:
    record(meta: Map, query: Path, subject: Path?, args: String?)

    output:
    record(meta: meta, tsv: file(subject ? 'query_vs_subject.blastn' : 'query.self.blastn'))

    script:
    def out = subject ? 'query_vs_subject.blastn' : 'query.self.blastn'
    """
    makeblastdb -in '${subject ?: query}' -dbtype nucl -out blastdb/subject
    blastn ${args ?: ''} -num_threads ${task.cpus} -query '${query}' -db blastdb/subject -outfmt "6 std qlen slen" -out '${out}'
    """

    stub:
    def out = subject ? 'query_vs_subject.blastn' : 'query.self.blastn'
    """
    touch '${out}'
    """
}
