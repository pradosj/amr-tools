nextflow.enable.types = true

process NCBI_DATASET_DOWNLOAD_GENOME {
    container 'docker.io/staphb/ncbi-datasets:18.18.0'
    memory 8.GB
    cpus 1
    time 30.m

    input:
    record(args: String)

    output:
    file('dataset')

    script:
    """
    datasets download genome ${args} --filename dataset.zip \\
        && unzip dataset.zip -d dataset
    """

    stub:
    """
    mkdir -p dataset/ncbi_dataset/data/GCF_000000000.1
    touch dataset/ncbi_dataset/data/assembly_data_report.jsonl
    touch dataset/ncbi_dataset/data/GCF_000000000.1/sequence_report.jsonl
    """
}
