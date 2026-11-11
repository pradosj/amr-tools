nextflow.enable.types = true

// The genes of a genome (FASTA file) to add to the pyMLST database, and its strain name
record IdFasta {
    path: Path
    id: String
}

// Whole/core genome MLST with pyMLST: create a database from the genes of `template`, add each genome (strain named by
// its id), then extract the distances, the MLST table and the multiple sequence alignment. `args` are given to
// `wgMLST add` (e.g. --coverage 0.85 --identity 0.85).
process PYMLST_WGMLST {
    container 'docker.io/staphb/pymlst:2.2.2'
    memory 8.GB
    cpus 1
    time 2.h

    input:
    record(meta: Record, template: Path, genomes: List<IdFasta>, args: String?)

    stage:
    stageAs(template, 'template.fasta')
    stageAs(genomes.collect({ g -> g.path }), 'genomes/?/*')    // the files may have the same name

    output:
    record(
        meta: meta,
        db: file('pymlst.db'),
        dir: file('pymlst')    // distance.tsv, mlst.tsv, msa.fasta
    )

    script:
    def adds = genomes.collect({ g -> "wgMLST add ${args ?: ''} --strain '${g.id}' pymlst.db '${g.path}'" }).join('\n')
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
