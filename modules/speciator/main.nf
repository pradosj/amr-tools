nextflow.enable.types = true

process SPECIATOR {
    container 'speciator-v4.0.0.sif'
    memory 10.GB
    cpus 1
    time 15.m

    input:
    record(meta: Map, fasta: Path)

    output:
    record(
        meta: meta,
        dir: file('speciator')
    )

    script:
    """
    mkdir ./speciator
    python3 /speciator.py ${fasta} /libraries /bactinspector/data/taxon_info.pqt > speciator/results.json
    """

    stub:
    """
    mkdir -p speciator
    touch speciator/results.json
    """
}
