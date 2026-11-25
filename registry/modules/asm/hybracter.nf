nextflow.enable.types = true

process HYBRACTER {
    container "quay.io/gbouras13/hybracter:0.11.2"
    memory 20.GB
    cpus 8
    time 4.h

    input:
    // long reads, plus optional short reads (short_read1, short_read2) for a hybrid assembly
    record(meta: Record?, short_read1: Path?, short_read2: Path?, long_reads: Path?, args: String?)

    output:
    // path: the final assembly, of FINAL_OUTPUT/complete (a contig longer than the chromosome length) or
    // FINAL_OUTPUT/incomplete; info: what the assembler tells of the contigs, by block: topology (circular or linear,
    // from the circular column of sample_per_contig_stats.tsv: only for the complete assemblies)
    record(meta: meta, path: file('assembly.fasta'), dir: file('hybracter'), info: file('info.json'))

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

    # the files of the sample (named sample by default) are in FINAL_OUTPUT/complete or FINAL_OUTPUT/incomplete
    final=hybracter/FINAL_OUTPUT/complete
    [ -f \$final/sample_final.fasta ] || final=hybracter/FINAL_OUTPUT/incomplete
    ln -s \$final/sample_final.fasta assembly.fasta

    # topology of each contig: circular column (True or False) of the per contig stats, missing for the incomplete assemblies
    awk -F'\\t' '
        BEGIN { printf "{\\n  \\"topology\\": {" }
        NR == 1 { for (i = 1; i <= NF; i++) if (\$i == "circular") c = i; next }
        c { printf "%s\\n    \\"%s\\": \\"%s\\"", (n++ ? "," : ""), \$1, (\$c == "True" ? "circular" : "linear") }
        END { print (n ? "\\n  " : "") "}\\n}" }
    ' \$final/sample_per_contig_stats.tsv > info.json
    """

    stub:
    """
    mkdir -p hybracter/FINAL_OUTPUT/complete
    touch hybracter/FINAL_OUTPUT/complete/sample_final.fasta
    ln -s hybracter/FINAL_OUTPUT/complete/sample_final.fasta assembly.fasta
    echo '{"topology": {}}' > info.json
    """
}
