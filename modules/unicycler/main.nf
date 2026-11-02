nextflow.enable.types = true

process UNICYCLER {
    container 'quay.io/biocontainers/unicycler:0.5.1--py312hdcc493e_4'
    errorStrategy task.attempt < 2 ? 'retry' : 'ignore'    // one retry with more resources, then ignore the sample
    memory 20.GB + 20.GB * (task.attempt - 1)
    cpus 8
    time 4.h + 8.h * (task.attempt - 1)

    input:
    // short reads (reads, read2) and/or long reads
    record(meta: Map, reads: Path?, read2: Path?, long_reads: Path?, args: String?)

    output:
    record(
        meta: meta,
        dir: file('unicycler'),
        fasta: file('unicycler/assembly.fasta')
    )

    script:
    def short_opts = !reads ? '' : (read2 ? "-1 ${reads} -2 ${read2}" : "-s ${reads}")
    def long_opts = long_reads ? "-l ${long_reads}" : ''
    """
    unicycler \\
        ${args ?: ''} \\
        --threads ${task.cpus} \\
        ${short_opts} ${long_opts} \\
        --out ./unicycler/
    """

    stub:
    """
    mkdir -p unicycler
    touch unicycler/assembly.fasta
    """
}
