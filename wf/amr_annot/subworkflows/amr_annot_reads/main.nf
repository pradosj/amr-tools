nextflow.enable.types = true

include { RESFINDER     as RESFINDER_LONGREAD      } from './modules/cgetools/resfinder'
include { RESFINDER     as RESFINDER_SHORTREAD     } from './modules/cgetools/resfinder'
include { PLASMIDFINDER as PLASMIDFINDER_LONGREAD  } from './modules/cgetools/plasmidfinder'
include { PLASMIDFINDER as PLASMIDFINDER_SHORTREAD } from './modules/cgetools/plasmidfinder'


// Options of a tool that can also run on the reads (params.amr.<tool>)
record ReadsToolParams {
	skip_longread: Boolean?    // skip on the long reads
	skip_shortread: Boolean?   // skip on the short reads
	args: String?
}

// Options of AMR_ANNOT_READS (the reads part of params.amr)
record AmrReadsParams {
	plasmidfinder: ReadsToolParams?
	resfinder: ReadsToolParams?
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


workflow AMR_ANNOT_READS {
	take:
	opts: AmrReadsParams
	short_reads: Channel<ShortReads>
	long_reads: Channel<LongReads>

	main:
	def plasmidfinder_long_ch = channel.empty()
	if (!opts.plasmidfinder?.skip_longread) {
		plasmidfinder_long_ch = PLASMIDFINDER_LONGREAD(
			long_reads.map({ r -> record(meta: r.meta, long_reads: r.reads, args: opts.plasmidfinder?.args) })
		)
	}
	def resfinder_long_ch = channel.empty()
	if (!opts.resfinder?.skip_longread) {
		resfinder_long_ch = RESFINDER_LONGREAD(
			long_reads.map({ r -> record(meta: r.meta, long_reads: r.reads, args: opts.resfinder?.args) })
		)
	}
	def plasmidfinder_short_ch = channel.empty()
	if (!opts.plasmidfinder?.skip_shortread) {
		plasmidfinder_short_ch = PLASMIDFINDER_SHORTREAD(short_reads.map({ r -> r + record(args: opts.plasmidfinder?.args) }))
	}
	def resfinder_short_ch = channel.empty()
	if (!opts.resfinder?.skip_shortread) {
		resfinder_short_ch = RESFINDER_SHORTREAD(short_reads.map({ r -> r + record(args: opts.resfinder?.args) }))
	}

	emit:
	resfinder_long: Channel<Record>      = resfinder_long_ch        // records with fields meta, dir
	resfinder_short: Channel<Record>     = resfinder_short_ch       // records with fields meta, dir
	plasmidfinder_long: Channel<Record>  = plasmidfinder_long_ch    // records with fields meta, dir
	plasmidfinder_short: Channel<Record> = plasmidfinder_short_ch   // records with fields meta, dir
}
