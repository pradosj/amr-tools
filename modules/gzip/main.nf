nextflow.enable.types = true

// Decompress a gzip file: X.fasta.gz -> X.fasta
process GZIP_DECOMPRESS {
    container 'docker.io/ubuntu:24.04'
    label 'local'
    memory 2.GB
    cpus 1
    time 30.m

    input:
    record(meta: Record, gz: Path)

    output:
    record(meta: meta, file: file(gz.baseName))

    script:
    """
    zcat '${gz}' > '${gz.baseName}'
    """

    stub:
    """
    touch '${gz.baseName}'
    """
}
