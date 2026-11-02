nextflow.enable.types = true

process PLASMIDFINDER {
    label 'cgetools'
    container "registry.gitlab.unige.ch/amr-genomics/cgetools:v2.0"
    memory 4.GB
    cpus 1
    time 1.h

    input:
    // either an assembly (fasta), short reads (reads, read2) or long reads
    record(meta: Map, fasta: Path?, reads: Path?, read2: Path?, long_reads: Path?, args: String?)

    output:
    record(
        meta: meta,
        dir: file('plasmidfinder')
    )

    script:
    def input = fasta ? "${fasta}"
        : long_reads ? "${long_reads}"
        : read2 ? "${reads} ${read2}"
        : "${reads}"
    """
    mkdir -p plasmidfinder
    plasmidfinder.py ${args ?: ''} -q -x -p /db/plasmidfinder_db -i ${input} -o 'plasmidfinder/'
    """

    stub:
    """
    mkdir -p plasmidfinder
    """
}
