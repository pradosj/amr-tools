nextflow.enable.types = true

process CHECKM2_DB {
    container 'quay.io/biocontainers/checkm2:1.1.0--pyh7e72e81_1'
    memory 4.GB
    cpus 1
    time 1.h

    output:
    file('checkm2_db')

    script:
    """
    checkm2 database --download --no_write_json_db --path checkm2_db
    """

    stub:
    """
    mkdir -p checkm2_db/CheckM2_database
    touch checkm2_db/CheckM2_database/uniref100.KO.1.dmnd
    """
}
