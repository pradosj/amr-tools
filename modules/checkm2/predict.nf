nextflow.enable.types = true

process CHECKM2_PREDICT {
    container 'quay.io/biocontainers/checkm2:1.1.0--pyh7e72e81_1'
    errorStrategy task.attempt < 2 ? 'retry' : 'ignore'    // one retry with more resources, then ignore the sample
    memory 8.GB + 8.GB * (task.attempt - 1)
    cpus 4
    time 30.m + 1.h * (task.attempt - 1)

    input:
    // db: CheckM2 database folder, as created by CHECKM2_DB
    record(meta: Record, fasta: Path, db: Path, args: String?)

    output:
    record(
        meta: meta,
        report: file('output/quality_report.tsv')
    )

    script:
    """
    checkm2 predict \\
        --input '${fasta}' \\
        --output-directory 'output' \\
        --threads ${task.cpus} \\
        --database_path '${db}/CheckM2_database/uniref100.KO.1.dmnd' \\
        ${args ?: ''}
    """

    stub:
    """
    mkdir -p output/ && touch output/quality_report.tsv
    """
}
