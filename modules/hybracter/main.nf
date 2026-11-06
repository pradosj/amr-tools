nextflow.enable.types = true

process HYBRACTER {
    container "quay.io/gbouras13/hybracter:0.11.2"
    memory 20.GB
    cpus 8
    time 4.h

    input:
    // long reads, plus optional short reads (short_read1, short_read2) for a hybrid assembly
    record(meta: Record, short_read1: Path?, short_read2: Path?, long_reads: Path?, args: String?)

    output:
    record(
        meta: meta,
        dir: file('hybracter'),
        fasta: file('hybracter/FINAL_OUTPUT/complete/sample_final.fasta')
    )

    script:
    def cmd = short_read1 ? 'hybrid-single' : 'long-single'
    def short_opts = !short_read1 ? '' : (short_read2 ? "-1 ${short_read1} -2 ${short_read2}" : "-s ${short_read1}")
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
