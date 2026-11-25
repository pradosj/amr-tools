nextflow.enable.types = true

process UNTAR {
    container 'docker.io/ubuntu:24.04'
    memory 4.GB
    cpus 2
    time 4.h

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

    stub:    // tree of the archive, with empty files (tar detects the compression when listing a file)
    """
    mkdir -p "${dir_name}"
    tar -tf "${archive}" | while read -r p; do
        case "\$p" in
            */) mkdir -p "${dir_name}/\$p" ;;
            *) mkdir -p "${dir_name}/\$(dirname "\$p")" && touch "${dir_name}/\$p" ;;
        esac
    done
    """
}
