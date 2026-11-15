nextflow.enable.types = true

process MULTIQC {
    container 'quay.io/biocontainers/multiqc:1.28--pyhdfd78af_0'
    cpus 2
    memory 2.GB
    time 1.h

    input:
    record(meta: Object, db: Path, config: Path, args: String?)

    output:
    record(
        meta: meta,
        html: file('multiqc.html')
    )

    script:
    """
    multiqc ${args ?: ''} --force --config ${config} --filename multiqc.html ${db}/
    """

    stub:
    """
    touch multiqc.html
    """
}
