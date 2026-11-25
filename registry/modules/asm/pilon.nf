nextflow.enable.types = true

process PILON {
    container 'quay.io/biocontainers/pilon:1.24--hdfd78af_0'
    memory 10.GB
    cpus 6
    time 4.h

    input:
    record(meta: Record?, path: Path, bam: Path, bai: Path, args: String?)    // path: the assembly to polish, bam: the short reads aligned on it

    output:
    record(meta: meta, path: file('pilon/pilon.fasta'), dir: file('pilon'))

    script:
    """
    mkdir pilon \\
    && java -Xmx8000M -jar /usr/local/share/pilon*/pilon.jar \\
        --genome ${path} \\
        --bam ${bam} \\
        ${args ?: '--fix all --changes'} \\
        --threads ${task.cpus} \\
        --output pilon/pilon
    """

    stub:
    """
    mkdir -p pilon
    touch pilon/pilon.fasta
    """
}
