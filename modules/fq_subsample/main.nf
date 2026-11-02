
process FQ_SUBSAMPLE {
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    memory '2 GB'
    cpus 1
    time '30 min'
    input:
    		tuple val(meta), path(reads)
    output:
    		tuple val(meta), path("subsampled_reads_*.fastq.gz")
    script:
    		reads = reads instanceof List?reads:[reads]
    		def bp_limit = task.ext.bp_limit?:-1
    		if (bp_limit<0) {
    			reads.withIndex().collect({x,i -> "ln -s '${x}' 'subsampled_reads_${i+1}.fastq.gz'"}).join("\n")
    		} else {
					"""
					set +o pipefail
					fq_subsample.awk -v bp_limit=${bp_limit} ${reads.join(" ")}
					"""
    		}
    stub:
    	reads.withIndex().collect({_x,i -> "touch 'subsampled_reads_${i+1}.fastq.gz'"}).join("\n")
}

