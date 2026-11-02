nextflow.enable.types = true

process HYBRACTER {
    container "quay.io/gbouras13/hybracter:0.11.2"
    memory 20.GB
    cpus 8
    time 4.h

    input:
    // long reads, plus optional short reads (reads, read2) for a hybrid assembly
    record(meta: Map, reads: Path?, read2: Path?, long_reads: Path?, args: String?)

    output:
    record(
        meta: meta,
        dir: file('hybracter'),
        fasta: file('hybracter/FINAL_OUTPUT/complete/sample_final.fasta')
    )

    script:
    def cmd = reads ? 'hybrid-single' : 'long-single'
    def short_opts = !reads ? '' : (read2 ? "-1 ${reads} -2 ${read2}" : "-s ${reads}")
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
