nextflow.enable.types = true

process AMRFINDERPLUS_UPDATE {
    container 'quay.io/biocontainers/ncbi-amrfinderplus:4.0.22--hf69ffd2_0'
    memory 4.GB
    cpus 1
    time 1.h

    output:
    file('amrfinder_db')

    script:
    """
    amrfinder_update -d db && mv \$(realpath db/latest) amrfinder_db
    """

    stub:
    """
    mkdir -p amrfinder_db
    """
}
