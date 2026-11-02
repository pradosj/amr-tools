nextflow.enable.types = true

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

params {
	gtf: String
}

process GTF2GFF {
	container "registry.gitlab.unige.ch/amr-genomics/rscript:v3"
	label 'rscript'
	memory 4.GB
	cpus 1
	time 30.m

	input:
	record(meta: Map, gtf: Path, assets: Path)

	stage:
	stageAs(gtf, 'ref.gtf')
	stageAs(assets, 'assets')

	output:
	record(meta: meta, gff: file('ref.gff'))

	script:
	"""
	#!/usr/bin/env Rscript
	source("assets/lib_gtf2gff.R")
	rtracklayer::import.gff2("ref.gtf") %>%
		gtf2gff() %>%
		rtracklayer::export.gff3("ref.gff")
	"""

	stub:
	"""
	touch ref.gff
	"""
}


workflow {
	main:
	// Validate parameters and print summary of supplied ones
	validateParameters()
	log.info(paramsSummaryLog(workflow))

	def assets = file("${projectDir}/assets")
	def gtf_ch = channel.fromPath(params.gtf).map({ f ->
		def id = f.name.replaceAll(/\.gtf(\.gz)?$/, '')
		record(meta: [sample_id: id], gtf: f, assets: assets)
	})

	def gff_ch = GTF2GFF(gtf_ch)

	publish:
	gff = gff_ch
}

output {
	gff: Channel<Record> {
		path { r -> r.gff >> "${r.meta.sample_id}.gff" }
	}
}
