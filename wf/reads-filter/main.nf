nextflow.enable.types = true

include { READS_FILTER; ReadsFilterParams } from './workflows/reads-filter/main.nf'
include { INPUT_READS; ReadsParams                        } from './subworkflows/input/main.nf'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'


// Input reads given as glob patterns (the sample_id is taken from the file names), and READS_FILTER options
params {
	reads: ReadsParams = record(ont: null, illumina: null)
	reads_filter: ReadsFilterParams = record(long_limit: 1_000_000_000, short_limit: 1_000_000_000, chopper_args: '--quality 15 --minlength 200')
}

workflow {
	main:
	// Validate parameters and print summary of supplied ones
	validateParameters()
	log.info(paramsSummaryLog(workflow))

	def reads = INPUT_READS(params.reads?.ont, params.reads?.illumina)
	def filtered = READS_FILTER(params.reads_filter, reads.short_reads, reads.long_reads)

	publish:
	reads_long_filtered = filtered.long_filtered
	reads_short_filtered = filtered.short_filtered
}

output {
	reads_long_filtered: Channel<Record> {
		path { r -> r.reads >> "reads/filtered_long/${r.meta.sample_id}.fastq.gz" }
	}
	reads_short_filtered: Channel<Record> {
		path { r ->
			r.reads >> "reads/filtered_short/${r.meta.sample_id}_R1.fastq.gz"
			r.read2 >> "reads/filtered_short/${r.meta.sample_id}_R2.fastq.gz"
		}
	}
}
