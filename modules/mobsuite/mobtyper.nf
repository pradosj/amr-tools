nextflow.enable.types = true

process MOBTYPER {
    container "quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_0"
    errorStrategy task.attempt < 2 ? 'retry' : 'ignore'    // one retry with more resources, then ignore the sample
    memory 16.GB + 16.GB * (task.attempt - 1)
    cpus 4
    time 30.m + 1.h * (task.attempt - 1)

    input:
    record(meta: Map, fasta: Path, args: String?)

    output:
    record(
        meta: meta,
        tsv: file('mobtyper.tsv')
    )

    script:
    """
    mob_typer ${args ?: ''} --num_threads=${task.cpus} --multi -i '${fasta}' -o mobtyper.tsv
    """

    stub:
    """
    touch mobtyper.tsv
    """
}
