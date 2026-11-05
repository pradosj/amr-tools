nextflow.enable.types = true

// Give each sequence of a FASTA file the id <prefix>_<n>, n being its rank in the file (headers reduced to the id),
// e.g. to keep the ids distinct when merging the sequences of several genomes
process FASTA_RENAME {
    container 'docker.io/ubuntu:24.04'
    label 'local'
    memory 2.GB
    cpus 1
    time 30.m

    input:
    record(meta: Map, fasta: Path, prefix: String)

    output:
    record(
        meta: meta,
        fasta: file('output.fasta'),    // the renamed sequences
        tsv:   file('map.tsv')          // original_id, new_id of each sequence
    )

    script:
    """
    awk -v prefix='${prefix}' '
    BEGIN { print "original_id\\tnew_id" > "map.tsv" }
    /^>/ {
        n++
        original_id = substr(\$0, 2)
        sub(/[ \\t].*\$/, "", original_id)
        new_id = prefix "_" n
        print ">" new_id
        print original_id "\\t" new_id > "map.tsv"
        next
    }
    { print }
    ' '${fasta}' > output.fasta
    """

    stub:
    """
    touch output.fasta map.tsv
    """
}
