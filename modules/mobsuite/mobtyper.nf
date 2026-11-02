nextflow.enable.types = true

process MOBTYPER {
    container "quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_0"
    memory 10.GB
    cpus 4
    time 30.m

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
