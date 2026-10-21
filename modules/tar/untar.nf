process UNTAR {
    container 'docker.io/ubuntu:24.04'
    memory '4 GB'
    cpus 2
    time '1h'
    input:
      tuple(path(archive),val(dirName))
    output:
      path "${dirName}/"
    script:
			def args = archive.name.endsWith('.tar') ? '' :
        (archive.name.endsWith('.tar.gz') || archive.name.endsWith('.tgz')) ? 'z' :
        (archive.name.endsWith('.tar.bz2')) ? 'j' :
        (archive.name.endsWith('.tar.xz')) ? 'J' :
        ''
        
	    """
	    mkdir -p "${dirName}"
	    tar -x${args} -f "${archive}" -C "${dirName}"
	    """
}
