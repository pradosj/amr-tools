nextflow.enable.types = true

include { FQ_SUBSAMPLE as FQ_SUBSAMPLE_LONG  } from './modules/fq_subsample'
include { CHOPPER as CHOPPER_LONG            } from './modules/chopper'
include { FQ_SUBSAMPLE as FQ_SUBSAMPLE_SHORT } from './modules/fq_subsample'


workflow READS_FILTER {
	take:
	opts: Map                 // params.reads_filter: [long_limit, short_limit, chopper_args]
	fqs_ch: Channel<Tuple>    // short reads: [meta, [R1, R2?]]
	fql_ch: Channel<Tuple>    // long reads:  [meta, fastq]

	main:
	def chopper_ch = CHOPPER_LONG(fql_ch, opts.chopper_args)
	def long_ch = FQ_SUBSAMPLE_LONG(
		chopper_ch.map({ r -> tuple(r.meta, [r.reads]) }),
		opts.long_limit as Long
	)
	def short_ch = FQ_SUBSAMPLE_SHORT(fqs_ch, opts.short_limit as Long)

	emit:
	long_filtered: Channel<Tuple>  = long_ch.map({ r -> tuple(r.meta, r.reads[0]) })
	short_filtered: Channel<Tuple> = short_ch.map({ r -> tuple(r.meta, r.reads) })
}
