nextflow.enable.types = true

process FASTP {
    container 'community.wave.seqera.io/library/fastp:1.0.1--c8b87fe62dcc103c'
    memory 4.GB
    cpus 4
    time 30.m

    input:
    record(meta: Record, read1: Path, read2: Path?, args: String?)    // ShortReads, plus args

    output:
    record(
        meta: meta,
        json: file('output.fastp.json'),
        html: file('output.fastp.html')
    )

    script:
    def in2 = read2 ? "--in2 ${read2}" : ''
    """
    fastp \\
        --in1 ${read1} ${in2} \\
        --thread ${task.cpus} \\
        --json output.fastp.json \\
        --html output.fastp.html \\
        ${args ?: ''}
    """

    stub:
    """
    touch output.fastp.json output.fastp.html
    """
}
