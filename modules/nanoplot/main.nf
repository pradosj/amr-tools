nextflow.enable.types = true

process NANOPLOT {
    container 'quay.io/biocontainers/nanoplot:1.41.0--pyhdfd78af_0'
    cpus 2
    memory 4.GB
    time 1.h

    input:
    tuple(meta: Map, reads: Path)
    args: String?

    output:
    record(
        meta: meta,
        nanoplot: file('nanoplot'),
        nanostat: file('nanoplot/NanoStats.txt')
    )

    script:
    """
    NanoPlot \\
        --threads ${task.cpus} \\
        --fastq ${reads} \\
        -o nanoplot \\
        ${args ?: ''}
    """

    stub:
    """
    mkdir -p nanoplot
    touch nanoplot/NanoStats.txt
    """
}
