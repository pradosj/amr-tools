nextflow.enable.types = true

// Sequences to align, and their meta: fasta is a FASTA file, or (subject only) the folder of a BLAST database made by
// MAKEBLASTDB (<folder>/db.*)
record BlastSequences {
    meta: Record
    fasta: Path
}

// BLASTN of the query sequences against the subject sequences, or against themselves when subject is null (do not give
// the same file as query and subject: it would be staged only once).
// A subject FASTA is made into a BLAST database to run on several threads (-num_threads is ignored with -subject).
// Output: tabular format with the query and subject lengths (-outfmt "6 std qlen slen"),
// query.self.blastn (self alignment) or query_vs_subject.blastn (empty when the query is empty).
process BLASTN {
    container 'docker.io/ncbi/blast:2.13.0'
    memory 8.GB
    cpus 4
    time 2.h

    input:
    record(query: BlastSequences, subject: BlastSequences?, args: String?)

    stage:
    stageAs(query.fasta, 'query.fasta')
    stageAs(subject?.fasta, 'subject')

    output:
    record(
        meta_query: query.meta,
        meta_subject: subject ? subject.meta : query.meta,    // the query meta for a self alignment
        tsv: file(subject ? 'query_vs_subject.blastn' : 'query.self.blastn')
    )

    script:
    def out = subject ? 'query_vs_subject.blastn' : 'query.self.blastn'
    """
    db=blastdb/db
    if [ -d subject ]; then
        db=subject/db
    else
        makeblastdb -in ${subject ? 'subject' : 'query.fasta'} -dbtype nucl -out "\$db"
    fi
    blastn ${args ?: ''} -num_threads ${task.cpus} -query query.fasta -db "\$db" -outfmt "6 std qlen slen" -out '${out}'
    """

    stub:
    def out = subject ? 'query_vs_subject.blastn' : 'query.self.blastn'
    """
    touch '${out}'
    """
}
