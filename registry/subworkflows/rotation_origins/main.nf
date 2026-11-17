nextflow.enable.types = true

include { FASTA_CAT             } from './modules/fasta/cat'
include { MAKEBLASTDB           } from './modules/blast/makeblastdb'
include { BLASTN as BLASTN_SELF } from './modules/blast/blastn'
include { IdPath; MetaPath      } from './subworkflows/input/main.nf'


// The regions shared by most of the genomes (at most one per sequence), from their self BLAST
process FIND_ROTATION_ORIGINS {
	container 'registry.gitlab.unige.ch/amr-genomics/rscript:v3'
	label 'rscript'
	memory 8.GB
	cpus 1
	time 30.m

	input:
	record(meta: Record?, fasta: Path, blast: Path)
	assets: Path    // the assets folder of the subworkflow (lib_io.R, lib_blast.R, lib_rotation_origins.R)

	stage:
	stageAs(fasta, 'genomes.fasta')
	stageAs(blast, 'genomes.self.blastn')
	stageAs(assets, 'assets')

	output:
	record(meta: meta, path: file('origins.fasta'))

	script:
	"""
	#!/usr/bin/env Rscript
	source("assets/lib_io.R")
	source("assets/lib_blast.R")
	source("assets/lib_rotation_origins.R")
	find_rotation_origins("genomes.fasta", "genomes.self.blastn", "origins.fasta")
	"""

	stub:
	"""
	touch origins.fasta
	"""
}


// A single FASTA of rotation origins for a set of genomes (one FASTA file per genome), to give as `origins` to
// GENOME_ROTATE: the longest regions shared by most of the genomes, at most one per sequence (chromosome, plasmid).
// Empty channel when there is no genome; empty FASTA when no region is shared.
workflow ROTATION_ORIGINS {
	take:
	meta: Record?                 // meta of the origins
	genomes: Channel<IdPath>      // the genomes (unique ids)
	blastn_args: String?          // self alignment of the genomes

	main:
	def fasta_ch = FASTA_CAT(genomes
		.collect()
		.flatMap({ rs -> rs.isEmpty() ? [] : [record(meta: meta, genomes: rs.toList().toSorted({ r -> r.meta.id }))] })
	).fasta
	def db_ch = MAKEBLASTDB(fasta_ch.map({ r -> record(meta: r.meta, path: r.path, dbtype: 'nucl', args: null) }))
	def self_ch = BLASTN_SELF(
		fasta_ch
			.combine(db_ch)
			.map({ f, d -> record(meta: f.meta, path: f.path, ref: d, args: blastn_args ?: '-perc_identity 70 -evalue 1e-3 -subject_besthit') })
	)
	def origins_ch = FIND_ROTATION_ORIGINS(
		fasta_ch
			.combine(self_ch)
			.map({ f, b -> record(meta: f.meta, fasta: f.path, blast: b.path) }),
		file("${moduleDir}/assets")
	)

	emit:
	origins: Channel<MetaPath> = origins_ch    // the FASTA of the origins (meta: the given meta), one record
}
