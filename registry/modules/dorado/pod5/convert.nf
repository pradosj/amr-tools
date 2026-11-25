
process DORADO_POD5_CONVERT {
    container 'docker.io/nanoporetech/dorado:shaf2aed69855de85e60b363c9be39558ef469ec365'
    memory '6 GB'
    cpus 4
    time '30 min'
    input:
    		tuple val(meta), path("input.fast5")
    output:
    		tuple val(meta), path("output.pod5")
    script:
    """
    pod5 convert from_fast5 -t ${task.cpus} input.fast5 --output output.pod5
    """
}


