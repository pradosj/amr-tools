nextflow.enable.types = true

// ANI of the query genomes against the reference genomes, or of the query genomes against each other when refs is empty
process FASTANI {
    container 'docker.io/staphb/fastani:1.34'
    memory 20.GB
    cpus 8
    time 15.m

    input:
    record(meta: Record, query: List<Path>, refs: List<Path>, args: String?)

    stage:
    stageAs(query, 'query/*')
    stageAs(refs, 'ref/*')

    output:
    record(meta: meta, tsv: file('fastani.tsv'))

    script:
    """
    find query/ -type f -or -type l > queryList.txt
    ${refs ? "find ref/ -type f -or -type l > refList.txt" : "cp queryList.txt refList.txt"}
    fastANI --threads ${task.cpus} ${args ?: ''} --refList refList.txt --queryList queryList.txt --output fastani.tsv
    """

    stub:
    """
    touch fastani.tsv
    """
}
