nextflow.enable.types = true

process AMRFINDERPLUS_RUN {
    container 'quay.io/biocontainers/ncbi-amrfinderplus:4.0.22--hf69ffd2_0'
    memory 6.GB
    cpus 4
    time 30.m

    input:
    // db: AMRFinderPlus database folder, as created by AMRFINDERPLUS_UPDATE
    record(meta: Record, fasta: Path, db: Path, args: String?)

    output:
    record(
        meta: meta,
        dir: file('amrfinderplus')
    )

    script:
    """
    mkdir amrfinderplus
    amrfinder ${args ?: ''} -n ${fasta} --mutation_all amrfinderplus/mutations.tsv --database ${db} --threads ${task.cpus} > amrfinderplus/report.tsv
    """

    stub:
    """
    mkdir -p amrfinderplus
    touch amrfinderplus/report.tsv amrfinderplus/mutations.tsv
    """
}
