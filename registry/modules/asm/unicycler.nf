nextflow.enable.types = true

process UNICYCLER {
    container 'quay.io/biocontainers/unicycler:0.5.1--py312hdcc493e_4'
    errorStrategy task.attempt < 2 ? 'retry' : 'ignore'    // one retry with more resources, then ignore the sample
    memory 20.GB + 20.GB * (task.attempt - 1)
    cpus 8
    time 4.h + 8.h * (task.attempt - 1)

    input:
    // short reads (short_read1, short_read2) and/or long reads
    record(meta: Record?, short_read1: Path?, short_read2: Path?, long_reads: Path?, args: String?)

    output:
    // asm_info: what the assembler tells of the contigs, by block: topology (circular or linear, from the FASTA headers)
    record(meta: meta, path: file('unicycler/assembly.fasta'), dir: file('unicycler'), asm_info: file('asm.json'))

    script:
    def short_opts = !short_read1 ? '' : (short_read2 ? "-1 ${short_read1} -2 ${short_read2}" : "-s ${short_read1}")
    def long_opts = long_reads ? "-l ${long_reads}" : ''
    """
    unicycler \\
        ${args ?: ''} \\
        --threads ${task.cpus} \\
        ${short_opts} ${long_opts} \\
        --out ./unicycler/

    # topology of each contig: circular=true in its FASTA header (e.g. >1 length=5329000 depth=1.00x circular=true), else linear
    awk '
        BEGIN { printf "{\\n  \\"topology\\": {" }
        /^>/ { printf "%s\\n    \\"%s\\": \\"%s\\"", (n++ ? "," : ""), substr(\$1, 2), (\$0 ~ / circular=true( |\$)/ ? "circular" : "linear") }
        END { print (n ? "\\n  " : "") "}\\n}" }
    ' unicycler/assembly.fasta > asm.json
    """

    stub:
    """
    mkdir -p unicycler
    touch unicycler/assembly.fasta
    echo '{"topology": {}}' > asm.json
    """
}
