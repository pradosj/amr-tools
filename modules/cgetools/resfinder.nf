nextflow.enable.types = true

process RESFINDER {
    label 'cgetools'
    container "registry.gitlab.unige.ch/amr-genomics/cgetools:v2.0"
    memory 4.GB
    cpus 3
    time 1.h

    input:
    // either an assembly (fasta), short reads (short_read1, short_read2) or long reads
    record(meta: Map, fasta: Path?, short_read1: Path?, short_read2: Path?, long_reads: Path?, args: String?)

    output:
    record(
        meta: meta,
        dir: file('resfinder')
    )

    script:
    def input = fasta ? "-ifa ${fasta}"
        : long_reads ? "--nanopore -ifq ${long_reads}"
        : short_read2 ? "-ifq ${short_read1} ${short_read2}"
        : "-ifq ${short_read1}"
    """
    python -m resfinder ${args ?: ''} ${input} --kma_threads ${task.cpus} -acq -d -j resfinder/data.json -o 'resfinder/'
    """

    stub:
    """
    mkdir -p resfinder
    """
}
