nextflow.enable.types = true

process MLST {
    label 'cgetools'
    container "registry.gitlab.unige.ch/amr-genomics/cgetools:main"
    memory 4.GB
    cpus 1
    time 1.h

    input:
    // species_args: MLST scheme selection, e.g. '--species ecoli'
    record(meta: Record, fasta: Path, species_args: String, args: String?)

    output:
    record(
        meta: meta,
        dir: file('mlst')
    )

    script:
    """
    mkdir -p mlst/ tmp/
    mlst.py -p ${args ?: ''} '/db/mlst_db' -i '${fasta}' -o 'mlst/' ${species_args} --tmp_dir 'tmp/' -x
    """

    stub:
    """
    mkdir -p mlst
    """
}
