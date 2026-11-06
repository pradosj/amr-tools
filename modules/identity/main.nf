nextflow.enable.types = true

// Copy input files into the work directory as task outputs (only task outputs can be published)
process IDENTITY {
    memory 2.GB
    cpus 1
    time 1.h
    label 'local'

    input:
    record(meta: Record, file1: Path, file2: Path?)

    output:
    record(
        meta: meta,
        file1: file(file1.name),
        file2: file(file2 ? file2.name : '.no_file2', optional: true)
    )

    script:
    "true"
}
