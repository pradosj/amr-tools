nextflow.enable.types = true

process MEDAKA_CONSENSUS {
    // the arm64 image (same tag) is set for label 'medaka' in the arm64 profile of conf/profiles.config
    container 'docker.io/ontresearch/medaka:shac4e11bfa4e65668b28739ba32edc3af12baf7574-amd64'
    label 'medaka'
    memory 10.GB
    cpus 4
    time 2.h

    input:
    record(meta: Record, path: Path, reads: Path, args: String?)    // path: the assembly to polish

    output:
    record(meta: meta, path: file('medaka/consensus.fasta'), dir: file('medaka'))

    script:
    """
    medaka_consensus \\
        ${args ?: '--bacteria'} \\
        -f -t ${task.cpus} \\
        -d ${path} \\
        -i ${reads} \\
        -o medaka
    """

    stub:
    """
    mkdir -p medaka/
    touch medaka/consensus.fasta
    """
}
