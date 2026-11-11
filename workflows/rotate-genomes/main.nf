nextflow.enable.types = true

include { MAKEBLASTDB } from './modules/blast/makeblastdb'
include { BLASTN as BLASTN_REFERENCE } from './modules/blast/blastn'

record MetaPath {
    meta: Record?
    path: Path
}

// The reference of a genome (same meta), one of: the reference dataset (FASTA file) to align on it, or the BLAST
// result of this alignment
record RotationReference {
	meta: Record
	fasta: Path?        // FASTA reference dataset, aligned with BLASTN on the genome
	blast: Path?        // hits of the reference sequences (query) on the genome (subject), -outfmt "6 std qlen slen"
}


// Rotate the sequences to start at their best hit of the reference, on its strand
process ROTATE_FASTA {
	container 'registry.gitlab.unige.ch/amr-genomics/rscript:v3'
	label 'rscript'
	memory 8.GB
	cpus 1
	time 30.m

	input:
	record(fasta: MetaPath, blast6: MetaPath)
	assets: Path    // the assets folder of the workflow (lib_io.R, lib_rotate_genomes.R)

	stage:
	stageAs(fasta.path, 'genome.fasta')
	stageAs(blast6.path, 'reference_vs_genome.blastn')
	stageAs(assets, 'assets')

	output:
	fasta: MetaPath = record(meta: fasta.meta, path: file('rotated.fasta'))

	script:
	"""
	#!/usr/bin/env Rscript
	source("assets/lib_io.R")
	source("assets/lib_rotate_genomes.R")
	rotate_fasta("genome.fasta", "reference_vs_genome.blastn", "rotated.fasta")
	"""

	stub:
	"""
	touch rotated.fasta
	"""
}


// Rotate genomes on a reference dataset: each sequence of a genome starts at its best hit of the reference,
// on the strand of the hit (sequences without hit are kept as they are). The hits come from the given BLAST
// result, or from the BLASTN of the reference on the genome. One output per genome: a genome without reference
// (same meta) or without hit is returned unchanged.
workflow ROTATE_GENOMES {
	take:
	fasta: Channel<MetaPath>
	references: Channel<RotationReference>
	blastn_args: String?    // BLASTN of the reference on the FASTA

	main:
	// Each FASTA paired with its reference: record(meta, fasta, ref), ref null when there is none (references without
	// FASTA are ignored)
	def input_ch = fasta
		.map({ f -> record(meta: f.meta, fasta: f) })
		.join(references.map({ r -> record(meta: r.meta, ref: r) }), by: 'meta', remainder: true)
		.filter({ r -> r.fasta != null })

	// Hits of each FASTA, from the first given of: the BLAST result, the FASTA reference
	def ref_fasta_ch = input_ch
		.filter({ r -> r.ref?.blast == null && r.ref?.fasta != null })
		.map({ r -> record(meta: r.meta, path: r.ref.fasta) })

	def genome_db_ch = MAKEBLASTDB(
		input_ch
			.join(ref_fasta_ch, by: 'meta')
			.map({ r -> record(fasta: r.fasta, args: null) })
	)
	def blast_ch = BLASTN_REFERENCE(
		ref_fasta_ch
			.join(genome_db_ch.map({ d -> record(meta: d.meta, db: d) }), by: 'meta')
			.map({ r -> record(query: record(meta: r.meta, path: r.path), subject_db: r.db, args: blastn_args) })
	)
	def hits_ch = input_ch
		.filter({ r -> r.ref?.blast != null })
		.map({ r -> record(meta: r.meta, path: r.ref.blast) })
		.mix(blast_ch)

	// One output per FASTA: the rotated FASTA, or the input FASTA when it has no reference or no hit
	def genome_hits_ch = input_ch
		.join(hits_ch.map({ r -> record(meta: r.meta, hits: r) }), by: 'meta', remainder: true)
	def rotated_ch = ROTATE_FASTA(
		genome_hits_ch
			.filter({ r -> r.hits != null && !(r.hits.path as Path).isEmpty() })
			.map({ r -> record(fasta: r.fasta, blast6: r.hits) }),
		file("${moduleDir}/assets")
	)
		.mix(genome_hits_ch
			.filter({ r -> r.hits == null || (r.hits.path as Path).isEmpty() })
			.map({ r -> r.fasta }))

	emit:
	rotated_fasta: Channel<MetaPath> = rotated_ch    // the rotated FASTA files (same meta)
	hits: Channel<MetaPath>          = hits_ch    // the hits of each genome (path: -outfmt "6 std qlen slen")
}
