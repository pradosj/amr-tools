nextflow.enable.types = true

// Archive the content of a folder (following the symbolic links); the compression is given by the
// extension of archive_name (.tar, .tar.gz/.tgz, .tar.bz2, .tar.xz)
process TAR {
    container 'docker.io/ubuntu:24.04'
    memory 4.GB
    cpus 2
    time 4.h

    input:
    record(dir: Path, archive_name: String)

    output:
    file(archive_name)

    script:
    def flag = (archive_name.endsWith('.tar.gz') || archive_name.endsWith('.tgz')) ? 'z' :
        archive_name.endsWith('.tar.bz2') ? 'j' :
        archive_name.endsWith('.tar.xz') ? 'J' :
        ''
    """
    tar -c${flag}h -f "${archive_name}" -C "${dir}" \$(ls -A "${dir}")
    """

    stub:
    """
    touch "${archive_name}"
    """
}
