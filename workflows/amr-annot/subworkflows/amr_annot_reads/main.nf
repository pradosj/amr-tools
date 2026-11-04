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
	meta: Map?
	read1: Path
	read2: Path?
}

// Long reads sample: one FASTQ per sample
record LongReads {
	meta: Map?
	reads: Path
}


// Gather the AMR results of the reads of a sample into a folder (published as amr/reads/<sample_id>.reads_amr), read into an RDS file
process AMR_COMBINE_READS {
	container 'registry.gitlab.unige.ch/amr-genomics/rscript:main'
	label 'rscript'
	memory 8.GB
	cpus 2
	time 30.m
	input:
	record(meta: Map, resfinder_long: Path?, resfinder_short: Path?, plasmidfinder_long: Path?, plasmidfinder_short: Path?)
	lib: Path    // lib_typing.R
	stage:
	stageAs(resfinder_long, 'reads_amr/resfinder_long')
	stageAs(resfinder_short, 'reads_amr/resfinder_short')
	stageAs(plasmidfinder_long, 'reads_amr/plasmidfinder_long')
	stageAs(plasmidfinder_short, 'reads_amr/plasmidfinder_short')
	output:
	record(meta: meta, dir: file('reads_amr'), rds: file('reads_amr.rds'))
	script:
	"""
	#!/usr/bin/env Rscript
	source("${lib}")
	read_reads_amr("reads_amr",'${meta.sample_id}') |> saveRDS("reads_amr.rds")
	"""
	stub:
	"""
	touch reads_amr.rds
	"""
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
		plasmidfinder_short_ch = PLASMIDFINDER_SHORTREAD(short_reads.map({ r -> record(meta: r.meta, short_read1: r.read1, short_read2: r.read2, args: opts.plasmidfinder?.args) }))
	}
	def resfinder_short_ch = channel.empty()
	if (!opts.resfinder?.skip_shortread) {
		resfinder_short_ch = RESFINDER_SHORTREAD(short_reads.map({ r -> record(meta: r.meta, short_read1: r.read1, short_read2: r.read2, args: opts.resfinder?.args) }))
	}

	// Gather the results of each sample into a folder and an RDS file (samples with at least one result)
	def combined_ch = AMR_COMBINE_READS(
		resfinder_long_ch.map({ r -> record(meta: r.meta, resfinder_long: r.dir) })
			.join(resfinder_short_ch.map({ r -> record(meta: r.meta, resfinder_short: r.dir) }), by: 'meta', remainder: true)
			.join(plasmidfinder_long_ch.map({ r -> record(meta: r.meta, plasmidfinder_long: r.dir) }), by: 'meta', remainder: true)
			.join(plasmidfinder_short_ch.map({ r -> record(meta: r.meta, plasmidfinder_short: r.dir) }), by: 'meta', remainder: true),
		file("${moduleDir}/assets/lib_typing.R")
	)

	emit:
	combined: Channel<Record>            = combined_ch              // records with fields meta, dir (results folder), rds (read_reads_amr() table)
	resfinder_long: Channel<Record>      = resfinder_long_ch        // records with fields meta, dir
	resfinder_short: Channel<Record>     = resfinder_short_ch       // records with fields meta, dir
	plasmidfinder_long: Channel<Record>  = plasmidfinder_long_ch    // records with fields meta, dir
	plasmidfinder_short: Channel<Record> = plasmidfinder_short_ch   // records with fields meta, dir
}
