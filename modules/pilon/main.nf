nextflow.enable.types = true

process PILON {
    container 'quay.io/biocontainers/pilon:1.24--hdfd78af_0'
    memory 10.GB
    cpus 6
    time 4.h

    input:
    record(meta: Record, fasta: Path, bam: Path, bai: Path, args: String?)

    output:
    record(
        meta: meta,
        dir: file('pilon'),
        fasta: file('pilon/pilon.fasta')
    )

    script:
    """
    mkdir pilon \\
    && java -Xmx8000M -jar /usr/local/share/pilon*/pilon.jar \\
        --genome ${fasta} \\
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
