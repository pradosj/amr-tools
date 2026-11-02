nextflow.enable.types = true

process SPADES {
    container 'quay.io/biocontainers/spades:4.2.0--h8d6e82b_1'
    memory 20.GB
    cpus 8
    time 4.h

    input:
    // short reads (reads, read2) and/or long reads
    record(meta: Map, reads: Path?, read2: Path?, long_reads: Path?, args: String?)

    output:
    record(
        meta: meta,
        dir: file('spades'),
        fasta: file('assembly.fasta')    // scaffolds, or contigs when there are no scaffolds
    )

    script:
    def short_opts = !reads ? '' : (read2 ? "-1 ${reads} -2 ${read2}" : "-s ${reads}")
    def long_opts = long_reads ? "--nanopore '${long_reads}'" : ''
    """
    mkdir -p spades && spades.py \\
        ${args ?: ''} \\
        --threads ${task.cpus} \\
        --memory ${task.memory.toGiga()} \\
        ${short_opts} ${long_opts} \\
        -o ./spades/
    if [ -f spades/scaffolds.fasta ]; then
        ln -s spades/scaffolds.fasta assembly.fasta
    else
        ln -s spades/contigs.fasta assembly.fasta
    fi
    """

    stub:
    """
    mkdir -p spades
    touch spades/scaffolds.fasta
    ln -s spades/scaffolds.fasta assembly.fasta
    """
}
