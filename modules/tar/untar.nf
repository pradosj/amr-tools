nextflow.enable.types = true

process UNTAR {
    container 'docker.io/ubuntu:24.04'
    memory 4.GB
    cpus 2
    time 1.h

    input:
    record(archive: Path, dir_name: String)

    output:
    file(dir_name)

    script:
    def name = archive.name
    def flag = (name.endsWith('.tar.gz') || name.endsWith('.tgz')) ? 'z' :
        name.endsWith('.tar.bz2') ? 'j' :
        name.endsWith('.tar.xz') ? 'J' :
        ''
    """
    mkdir -p "${dir_name}"
    tar -x${flag} -f "${archive}" -C "${dir_name}"
    """

    stub:
    """
    mkdir -p "${dir_name}"
    """
}
