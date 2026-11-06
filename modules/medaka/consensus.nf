nextflow.enable.types = true

process MEDAKA_CONSENSUS {
    // the arm64 image (same tag) is set for label 'medaka' in the arm64 profile of conf/profiles.config
    container 'docker.io/ontresearch/medaka:shac4e11bfa4e65668b28739ba32edc3af12baf7574-amd64'
    label 'medaka'
    memory 10.GB
    cpus 4
    time 2.h

    input:
    record(meta: Record, fasta: Path, reads: Path, args: String?)

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
