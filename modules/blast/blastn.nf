nextflow.enable.types = true

// BLASTN of the query sequences against the subject sequences, or against themselves when subject is null.
// Output: tabular format with the query and subject lengths (-outfmt "6 std qlen slen"), named after the
// query and subject files (<query>.self.blastn or <query>_vs_<subject>.blastn); empty when the query is empty.
process BLASTN {
    container 'docker.io/ncbi/blast:2.13.0'
    memory 8.GB
    cpus 1
    time 2.h

    input:
    record(meta: Map, query: Path, subject: Path?, args: String?)

    output:
    record(meta: meta, tsv: file(subject ? "${query.baseName}_vs_${subject.baseName}.blastn" : "${query.baseName}.self.blastn"))

    script:
    def out = subject ? "${query.baseName}_vs_${subject.baseName}.blastn" : "${query.baseName}.self.blastn"
    """
    if [ -s '${query}' ]; then
        blastn ${args ?: ''} -query '${query}' -subject '${subject ?: query}' -outfmt "6 std qlen slen" -out '${out}'
    else
        touch '${out}'
    fi
    """

    stub:
    def out = subject ? "${query.baseName}_vs_${subject.baseName}.blastn" : "${query.baseName}.self.blastn"
    """
    touch '${out}'
    """
}
