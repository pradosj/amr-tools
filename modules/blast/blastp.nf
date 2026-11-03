nextflow.enable.types = true

// BLASTP of the query sequences against the subject sequences, or against themselves when subject is null.
// Output: tabular format with the query and subject lengths (-outfmt "6 std qlen slen"), named after the
// query and subject files (<query>.self.blastp or <query>_vs_<subject>.blastp); empty when the query is empty.
process BLASTP {
    container 'docker.io/ncbi/blast:2.13.0'
    memory 8.GB
    cpus 1
    time 2.h

    input:
    record(meta: Map, query: Path, subject: Path?, args: String?)

    output:
    record(meta: meta, tsv: file(subject ? "${query.baseName}_vs_${subject.baseName}.blastp" : "${query.baseName}.self.blastp"))

    script:
    def out = subject ? "${query.baseName}_vs_${subject.baseName}.blastp" : "${query.baseName}.self.blastp"
    """
    if [ -s '${query}' ]; then
        blastp ${args ?: ''} -query '${query}' -subject '${subject ?: query}' -outfmt "6 std qlen slen" -out '${out}'
    else
        touch '${out}'
    fi
    """

    stub:
    def out = subject ? "${query.baseName}_vs_${subject.baseName}.blastp" : "${query.baseName}.self.blastp"
    """
    touch '${out}'
    """
}
