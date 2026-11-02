nextflow.enable.types = true

process FQ_SUBSAMPLE {
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    memory 2.GB
    cpus 1
    time 30.m

    input:
    record(meta: Map, reads: List<Path>, bp_limit: Long?)   // bp_limit: maximum total read length in bp (default 1_000_000_000), negative to keep all reads

    output:
    record(
        meta: meta,
        reads: files('subsampled_reads_*.fastq.gz').toSorted()  // glob is unordered: sort to match input order (subsampled_reads_<input index>)
    )

    script:
    def limit = bp_limit != null ? bp_limit : 1_000_000_000L
    if (limit < 0) {
        reads.withIndex().collect({ x, i -> "ln -s '${x}' 'subsampled_reads_${i+1}.fastq.gz'" }).join('\n')
    } else {
        """
        set +o pipefail
        fq_subsample.awk -v bp_limit=${limit} ${reads.join(' ')}
        """
    }

    stub:
    reads.withIndex().collect({ _x, i -> "touch 'subsampled_reads_${i+1}.fastq.gz'" }).join('\n')
}
