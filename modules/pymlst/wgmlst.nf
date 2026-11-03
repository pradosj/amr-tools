nextflow.enable.types = true

// Whole/core genome MLST with pyMLST: create a database from the genes of `template`, add each genome
// (strain named after the file, without extension), then extract the distances, the MLST table and the
// multiple sequence alignment. `args` are given to `wgMLST add` (e.g. --coverage 0.85 --identity 0.85).
process PYMLST_WGMLST {
    container 'docker.io/staphb/pymlst:2.2.2'
    memory 8.GB
    cpus 1
    time 2.h

    input:
    record(meta: Map, template: Path, genomes: List<Path>, args: String?)

    stage:
    stageAs(template, 'template.fasta')
    stageAs(genomes, 'genomes/*')

    output:
    record(
        meta: meta,
        db: file('pymlst.db'),
        dir: file('pymlst')    // distance.tsv, mlst.tsv, msa.fasta
    )

    script:
    """
    wgMLST create --force pymlst.db template.fasta
    for f in genomes/*; do
        wgMLST add ${args ?: ''} --strain "\$(basename "\${f%.*}")" pymlst.db "\$f"
    done
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
