nextflow.enable.types = true

process MEDAKA_CONSENSUS {
    container 'docker.io/ontresearch/medaka:shac4e11bfa4e65668b28739ba32edc3af12baf7574-amd64'
    memory 10.GB
    cpus 4
    time 2.h

    input:
    record(meta: Map, fasta: Path, reads: Path, args: String?)

    output:
    record(
        meta: meta,
        dir: file('medaka'),
        fasta: file('medaka/consensus.fasta')
    )

    script:
    """
    medaka_consensus \\
        ${args ?: '--bacteria'} \\
        -f -t ${task.cpus} \\
        -d ${fasta} \\
        -i ${reads} \\
        -o medaka
    """

    stub:
    """
    mkdir -p medaka/
    touch medaka/consensus.fasta
    """
}
