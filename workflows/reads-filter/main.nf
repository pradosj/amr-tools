nextflow.enable.types = true

include { FQ_SUBSAMPLE as FQ_SUBSAMPLE_LONG  } from './modules/fq_subsample'
include { CHOPPER as CHOPPER_LONG            } from './modules/chopper'
include { FQ_SUBSAMPLE as FQ_SUBSAMPLE_SHORT } from './modules/fq_subsample'


// Options of READS_FILTER (params.reads_filter). Fields are nullable: modules fall back to their own defaults.
record ReadsFilterParams {
	long_limit: Integer?       // FQ_SUBSAMPLE bp limit for long reads, negative to keep all reads
	short_limit: Integer?      // FQ_SUBSAMPLE bp limit for short reads, negative to keep all reads
	chopper_args: String?      // CHOPPER arguments
}

// Illumina reads of a sample: read2 is null for single-end reads
record ShortReads {
	meta: Map?
	read1: Path
	read2: Path?
}

// ONT reads of a sample: one FASTQ
record LongReads {
	meta: Map?
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
		chopper_ch.map({ r -> record(meta: r.meta, read1: r.reads, bp_limit: opts.long_limit as Long) })
	).map({ r -> record(meta: r.meta, reads: r.read1) })
	def short_ch = FQ_SUBSAMPLE_SHORT(
		short_reads.map({ r -> r + record(bp_limit: opts.short_limit as Long) })
	)

	emit:
	long_filtered: Channel<LongReads>   = long_ch
	short_filtered: Channel<ShortReads> = short_ch
}
