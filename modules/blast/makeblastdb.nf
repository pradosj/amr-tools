nextflow.enable.types = true

record MetaPath {
    meta: Record?
    path: Path
}

// BLAST database of the sequences of a FASTA file (dbtype: nucl, the default, or prot), for BLASTN/BLASTP `subject_db`
// Output: the folder of the database files, blastdb/db.*
process MAKEBLASTDB {
    container 'docker.io/ncbi/blast:2.13.0'
    memory 4.GB
    cpus 1
    time 30.m

    input:
    record(fasta: MetaPath, dbtype: String?, args: String?)    // dbtype: nucl or prot (nucl by default)

    output:
    db: MetaPath = record(meta: fasta.meta, path: file('blastdb'))

    script:
    """
    makeblastdb ${args ?: ''} -in '${fasta.path}' -dbtype ${dbtype ?: 'nucl'} -out blastdb/db
    """

    stub:
    """
    mkdir blastdb
    touch blastdb/db.pin
    """
}
