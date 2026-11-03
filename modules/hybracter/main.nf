nextflow.enable.types = true

process HYBRACTER {
    container "quay.io/gbouras13/hybracter:0.11.2"
    memory 20.GB
    cpus 8
    time 4.h

    input:
    // long reads, plus optional short reads (short_reads, short_reads2) for a hybrid assembly
    record(meta: Map, short_reads: Path?, short_reads2: Path?, long_reads: Path?, args: String?)

    output:
    record(
        meta: meta,
        dir: file('hybracter'),
        fasta: file('hybracter/FINAL_OUTPUT/complete/sample_final.fasta')
    )

    script:
    def cmd = short_reads ? 'hybrid-single' : 'long-single'
    def short_opts = !short_reads ? '' : (short_reads2 ? "-1 ${short_reads} -2 ${short_reads2}" : "-s ${short_reads}")
    def long_opts = long_reads ? "-l ${long_reads}" : ''
    """
    hybracter ${cmd} \\
        ${args ?: ''} \\
        --threads ${task.cpus} \\
        --auto \\
        ${short_opts} ${long_opts} \\
        -o hybracter
    """

    stub:
    """
    mkdir -p hybracter/FINAL_OUTPUT/complete
    touch hybracter/FINAL_OUTPUT/complete/sample_final.fasta
    """
}
