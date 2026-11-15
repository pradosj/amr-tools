nextflow.enable.types = true

// Download the CheckM2 database into the folder dir_name.
// dir_name is also needed as an input: a typed process without input and with a stub block is never cached
// (Nextflow 26.10.4 and 26.10.5), so the database would be downloaded again at each run.
process CHECKM2_DB {
    container 'quay.io/biocontainers/checkm2:1.1.0--pyh7e72e81_1'
    memory 4.GB
    cpus 1
    time 1.h

    input:
    dir_name: String

    output:
    file(dir_name)

    script:
    """
    checkm2 database --download --no_write_json_db --path '${dir_name}'
    """

    stub:
    """
    mkdir -p '${dir_name}/CheckM2_database'
    touch '${dir_name}/CheckM2_database/uniref100.KO.1.dmnd'
    """
}
