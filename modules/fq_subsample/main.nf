nextflow.enable.types = true

process FQ_SUBSAMPLE {
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    memory 2.GB
    cpus 1
    time 30.m

    input:
    record(meta: Map, reads: Path, read2: Path?, bp_limit: Long?)   // bp_limit: maximum total read length in bp (default 1_000_000_000), negative to keep all reads

    output:
    record(
        meta: meta,
        reads: file('subsampled_reads_1.fastq.gz'),
        read2: file('subsampled_reads_2.fastq.gz', optional: true)
    )

    script:
    def limit = bp_limit != null ? bp_limit : 1_000_000_000L
    def r2 = read2 ? "${read2}" : ''
    if (limit < 0) {
        """
        ln -s '${reads}' 'subsampled_reads_1.fastq.gz'
        ${read2 ? "ln -s '${read2}' 'subsampled_reads_2.fastq.gz'" : ''}
        """
    } else {
        """
        set +o pipefail
        fq_subsample.awk -v bp_limit=${limit} ${reads} ${r2}
        """
    }

    stub:
    """
    touch subsampled_reads_1.fastq.gz
    ${read2 ? 'touch subsampled_reads_2.fastq.gz' : ''}
    """
}
