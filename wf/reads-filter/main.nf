nextflow.enable.types = true

include { FQ_SUBSAMPLE as FQ_SUBSAMPLE_LONG  } from './modules/fq_subsample'
include { CHOPPER as CHOPPER_LONG            } from './modules/chopper'
include { FQ_SUBSAMPLE as FQ_SUBSAMPLE_SHORT } from './modules/fq_subsample'
include { READS_INPUT; ReadsParams                        } from './conf/inputs.nf'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'


// Options of READS_FILTER (params.reads_filter). Fields are nullable: modules fall back to their own defaults.
record ReadsFilterParams {
	long_limit: Integer?       // FQ_SUBSAMPLE bp limit for long reads, negative to keep all reads
	short_limit: Integer?      // FQ_SUBSAMPLE bp limit for short reads, negative to keep all reads
	chopper_args: String?      // CHOPPER arguments
}

// Short reads sample: read2 is null for single-end reads
record ShortReads {
	meta: Map
	reads: Path
	read2: Path?
}

// Long reads sample: one FASTQ per sample
record LongReads {
	meta: Map
	reads: Path
}

workflow READS_FILTER {
	take:
	opts: ReadsFilterParams
	short_reads: Channel<ShortReads>
	long_reads: Channel<LongReads>

	main:
	def chopper_ch = CHOPPER_LONG(
		long_reads.map({ r -> r + record(args: opts.chopper_args) })
	)
	def long_ch = FQ_SUBSAMPLE_LONG(
		chopper_ch.map({ r -> r + record(bp_limit: opts.long_limit as Long) })
	)
	def short_ch = FQ_SUBSAMPLE_SHORT(
		short_reads.map({ r -> r + record(bp_limit: opts.short_limit as Long) })
	)

	emit:
	long_filtered: Channel<LongReads>   = long_ch
	short_filtered: Channel<ShortReads> = short_ch
}


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

	def reads = READS_INPUT(params.reads?.ont, params.reads?.illumina)
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
