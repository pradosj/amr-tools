nextflow.enable.types = true

// BLAST database of the sequences of a FASTA file (dbtype: prot or nucl), for BLASTP/BLASTN `db`
// Output: the folder of the database files, blastdb/db.*
process MAKEBLASTDB {
    container 'docker.io/ncbi/blast:2.13.0'
    memory 4.GB
    cpus 1
    time 30.m

    input:
    record(meta: Record, fasta: Path, dbtype: String, args: String?)

    output:
    record(meta: meta, db: file('blastdb'))

    script:
    """
    makeblastdb ${args ?: ''} -in '${fasta}' -dbtype ${dbtype} -out blastdb/db
    """

    stub:
    """
    mkdir blastdb
    touch blastdb/db.pin
    """
}
