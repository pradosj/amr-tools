nextflow.enable.types = true

include { FASTA_CAT     } from './modules/fasta/cat'
include { MAKEBLASTDB   } from './modules/blast/makeblastdb'
include { BLASTN        } from './modules/blast/blastn'
include { QUARTO_RENDER } from './modules/quarto/render'
include { IdPath; MetaPath } from './subworkflows/input/main.nf'


process GENOME_DISTANCE {
	container 'registry.gitlab.unige.ch/amr-genomics/rscript:v3'
	label 'rscript'
	memory 8.GB
	cpus 1
	time 30.m

	input:
	record(meta: Record?, blast6: MetaPath, genomes_tsv: MetaPath, min_len: Integer)
	assets: Path

	stage:
	stageAs(blast6.path, 'genomes.self.blastn')
	stageAs(genomes_tsv.path, 'genomes.tsv')
	stageAs(assets, 'assets')

	output:
	record(meta: meta, graph: file('distance_graph.rds'))

	script:
	"""
	#!/usr/bin/env Rscript
	source("assets/lib_wg_blastn_distance.R")
	seq_genomes <- read_tsv("genomes.tsv", col_types = "ccc")
	g <- genome_distance_graph("genomes.self.blastn", seq_genomes, min_len = ${min_len})
	saveRDS(g, "distance_graph.rds")
	"""

	stub:
	"""
	touch distance_graph.rds
	"""
}


workflow WG_BLASTN_DISTANCE {
	take:
	genomes: Channel<IdPath>
	blastn_args: String?
	min_len: Integer?

	main:
	def cat_ch = FASTA_CAT(genomes.collect().map({ rs -> record(genomes: rs.toList().toSorted({ r -> r.meta.id })) }))
	def db_ch = MAKEBLASTDB(cat_ch.fasta)
	def blastn_ch = BLASTN(cat_ch.fasta.combine(ref: db_ch, args: blastn_args ?: '-perc_identity 70 -evalue 1e-3 -subject_besthit -max_target_seqs 1000000'))
	def distance_ch = GENOME_DISTANCE(
		channel.value(record(min_len: min_len ?: 1000)).combine(blast6: blastn_ch, genomes_tsv: cat_ch.tsv),
		file("${moduleDir}/assets")
	)

	def report_ch = QUARTO_RENDER(distance_ch.map({ r -> record(
		meta: r.meta,
		template: file("${moduleDir}/assets/wg_blastn_distance_report.qmd"),
		extra: null,
		files: [r.graph],
		args: "-P 'graph:distance_graph.rds'"
	) }))

	emit:
	distance = distance_ch
	report = report_ch
}
