nextflow.enable.types = true

// Coverage of each reference sequence (contig) of a BAM: numreads, covbases, coverage (%), meandepth, ...
process SAMTOOLS_COVERAGE {
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    memory 4.GB
    cpus 1
    time 30.m

    input:
    record(meta: Record?, path: Path, args: String?)

    output:
    coverage = record(
        meta: meta,
        path: file("${path.name}.coverage.tsv")
    )

    script:
    """
    samtools coverage ${args ?: ''} ${path} > ${path.name}.coverage.tsv
    """

    stub:
    """
    touch ${path.name}.coverage.tsv
    """
}
