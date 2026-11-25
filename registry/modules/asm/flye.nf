nextflow.enable.types = true

process FLYE {
    container 'quay.io/biocontainers/flye:2.9.6--py311h2de2dd3_0'
    errorStrategy task.attempt < 2 ? 'retry' : 'ignore'    // one retry with more resources, then ignore the sample
    memory 20.GB + 32.GB * (task.attempt - 1)
    cpus 8
    time 4.h + 8.h * (task.attempt - 1)

    input:
    record(meta: Record?, long_reads: Path, args: String?)

    output:
    // info: what the assembler tells of the contigs, by block: topology (circular or linear, from assembly_info.txt)
    record(meta: meta, path: file('flye/assembly.fasta'), dir: file('flye'), info: file('info.json'))

    script:
    """
    flye \\
        ${args ?: ''} \\
        --threads ${task.cpus} \\
        --out-dir ./flye \\
        --nano-hq ${long_reads}

    # topology of each contig: circ. column (4th) of assembly_info.txt, Y: circular, N: linear
    awk -F'\\t' '
        BEGIN { printf "{\\n  \\"topology\\": {" }
        NR > 1 { printf "%s\\n    \\"%s\\": \\"%s\\"", (n++ ? "," : ""), \$1, (\$4 == "Y" ? "circular" : "linear") }
        END { print (n ? "\\n  " : "") "}\\n}" }
    ' flye/assembly_info.txt > info.json
    """

    stub:
    """
    mkdir -p flye/
    touch flye/assembly.fasta
    echo '{"topology": {}}' > info.json
    """
}
