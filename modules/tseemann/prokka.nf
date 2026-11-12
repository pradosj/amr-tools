nextflow.enable.types = true

process PROKKA {
    container "docker.io/staphb/prokka:1.14.6"
    memory 12.GB
    cpus 4
    time 1.h

    input:
    record(meta: Record?, path: Path, args: String?)

    output:
    dir: record(meta: meta, path: file('prokka'))              // all the prokka outputs
    faa: record(meta: meta, path: file('prokka/prokka.faa'))   // proteins
    ffn: record(meta: meta, path: file('prokka/prokka.ffn'))   // nucleotide sequences of the genes
    gff: record(meta: meta, path: file('prokka/prokka.gff'))   // annotations, followed by the sequences
    fna: record(meta: meta, path: file('prokka/prokka.fna'))   // annotated sequences

    script:
    """
    prokka ${args ?: '--kingdom Bacteria'} --outdir 'prokka/' --cpus ${task.cpus} --prefix prokka ${path}
    """

    stub:
    """
    mkdir -p prokka
    touch prokka/prokka.faa prokka/prokka.ffn prokka/prokka.fna prokka/prokka.gff prokka/prokka.tsv
    """
}
