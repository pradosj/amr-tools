nextflow.enable.types = true

record FastaMeta {
	id: String
}

record IdFasta {
	meta: FastaMeta
	path: Path
}

process PYMLST_WGMLST {
    container 'docker.io/staphb/pymlst:2.2.2'
    memory 8.GB
    cpus 1
    time 2.h

    input:
    record(meta: Record?, template: Path, genomes: List<IdFasta>, args: String?)

    stage:
    stageAs(template, 'template.fasta')
    stageAs(genomes.collect({ g -> g.path }), 'genomes/?/*')

    output:
    record(
        meta: meta,
        db: file('pymlst.db'),
        dir: file('pymlst')
    )

    script:
    def adds = genomes.collect({ g -> "wgMLST add ${args ?: ''} --strain '${g.meta.id}' pymlst.db '${g.path}'" }).join('\n')
    """
    wgMLST create --force pymlst.db template.fasta
    ${adds}
    mkdir pymlst
    wgMLST distance --output pymlst/distance.tsv pymlst.db
    wgMLST mlst --output pymlst/mlst.tsv pymlst.db
    wgMLST msa --output pymlst/msa.fasta pymlst.db
    """

    stub:
    """
    touch pymlst.db
    mkdir pymlst
    touch pymlst/distance.tsv pymlst/mlst.tsv pymlst/msa.fasta
    """
}
