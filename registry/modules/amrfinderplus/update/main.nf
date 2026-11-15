nextflow.enable.types = true

// Download the latest AMRFinderPlus database into the folder dir_name.
// dir_name is also needed as an input: a typed process without input and with a stub block is never cached
// (Nextflow 26.10.4 and 26.10.5), so the database would be downloaded again at each run.
process AMRFINDERPLUS_UPDATE {
    container 'quay.io/biocontainers/ncbi-amrfinderplus:4.0.22--hf69ffd2_0'
    memory 4.GB
    cpus 1
    time 1.h

    input:
    dir_name: String

    output:
    file(dir_name)

    script:
    """
    amrfinder_update -d db && mv \$(realpath db/latest) '${dir_name}'
    """

    stub:
    """
    mkdir -p '${dir_name}'
    """
}
