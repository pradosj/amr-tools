nextflow.enable.types = true

// BLAST database of the sequences of a FASTA file (dbtype: nucl, the default, or prot), for BLASTN/BLASTP `ref`
// Output: the folder of the database files, blastdb/db.*
process MAKEBLASTDB {
    container 'docker.io/ncbi/blast:2.13.0'
    memory 4.GB
    cpus 1
    time 30.m

    input:
    record(meta: Record?, path: Path, dbtype: String?, args: String?)    // dbtype: nucl or prot (nucl by default)

    output:
    db = record(meta: meta, path: file('blastdb'))

    script:
    """
    makeblastdb ${args ?: ''} -in '${path}' -dbtype ${dbtype ?: 'nucl'} -out blastdb/db
    """

    stub:
    """
    mkdir blastdb
    touch blastdb/db.pin
    """
}
