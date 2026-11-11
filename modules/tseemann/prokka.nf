nextflow.enable.types = true

record MetaPath {
    meta: Record?
    path: Path
}

process PROKKA {
    container "docker.io/staphb/prokka:1.14.6"
    memory 12.GB
    cpus 4
    time 1.h

    input:
    record(fasta: MetaPath, args: String?)

    output:
    record(
        meta: fasta.meta,
        dir: file('prokka'),              // all the prokka outputs
        faa: file('prokka/prokka.faa'),   // proteins
        ffn: file('prokka/prokka.ffn'),   // nucleotide sequences of the genes
        gff: file('prokka/prokka.gff'),   // annotations, followed by the sequences
        fna: file('prokka/prokka.fna')    // annotated sequences
    )

    script:
    """
    prokka ${args ?: '--kingdom Bacteria'} --outdir 'prokka/' --cpus ${task.cpus} --prefix prokka ${fasta.path}
    """

    stub:
    """
    mkdir -p prokka
    touch prokka/prokka.faa prokka/prokka.ffn prokka/prokka.fna prokka/prokka.gff prokka/prokka.tsv
    """
}
