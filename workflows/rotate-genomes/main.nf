nextflow.enable.types = true

include { MAKEBLASTDB } from './modules/blast/makeblastdb'
include { BLASTN as BLASTN_REFERENCE } from './modules/blast/blastn'


// A genome: the sequences to rotate, or the rotated sequences (gff and info are not set on the rotated genomes)
record Genome {
	meta: Record
	fasta: Path
	gff: Path?
	info: Path?
}

// The reference of a genome (same meta), one of: the reference dataset (FASTA file or string) to align on it, or the
// BLAST result of this alignment
record RotationReference {
	meta: Record
	fasta: Path?        // FASTA reference dataset, aligned with BLASTN on the genome
	sequence: String?   // reference as a string: FASTA text (starting with ">"), or a sequence (written under the header ">reference")
	blast: Path?        // hits of the reference sequences (query) on the genome (subject), -outfmt "6 std qlen slen"
}

// Hits of the reference sequences (query) on a genome (subject), -outfmt "6 std qlen slen"
record ReferenceHits {
	meta: Record
	blast6: Path
}


// Reference given as a string: FASTA text (starting with '>'), or a sequence written under the header '>reference'
process WRITE_FASTA {
	input:
	record(meta: Record, sequence: String)
	output:
	record(meta: meta, fasta: file('reference.fasta'))
	exec:    // native task: the file is written by Nextflow, without a shell script
	def fasta = sequence.trim().startsWith('>') ? sequence.trim() : ">reference\n${sequence.replaceAll(/\s+/, '')}"
	task.workDir.resolve('reference.fasta').text = fasta + '\n'
}


// Rotate the sequences to start at their best hit of the reference, on its strand
process ROTATE_FASTA {
	container 'registry.gitlab.unige.ch/amr-genomics/rscript:v3'
	label 'rscript'
	memory 8.GB
	cpus 1
	time 30.m

	input:
	record(meta: Record, fasta: Path, blast: Path)
	assets: Path    // the assets folder of the workflow (lib_io.R, lib_rotate_genomes.R)

	stage:
	stageAs(fasta, 'genome.fasta')
	stageAs(blast, 'reference_vs_genome.blastn')
	stageAs(assets, 'assets')

	output:
	record(meta: meta, fasta: file('rotated.fasta'))

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
	genomes: Channel<Genome>
	references: Channel<RotationReference>
	blastn_args: String?    // BLASTN of the reference on the genomes

	main:
	// Each genome paired with its reference: record(meta, genome, ref), ref null when there is none
	// (references without genome are ignored)
	def input_ch = genomes
		.map({ g -> record(meta: g.meta, genome: g) })
		.join(references.map({ r -> record(meta: r.meta, ref: r) }), by: 'meta', remainder: true)
		.filter({ r -> r.genome != null })

	// Hits of each genome, from the first given of: the BLAST result, the FASTA reference, the reference sequence
	def ref_fasta_ch = input_ch
		.filter({ r -> r.ref?.blast == null && r.ref?.fasta != null })
		.map({ r -> record(meta: r.meta, ref_fasta: r.ref.fasta) })
		.mix(WRITE_FASTA(
			input_ch
				.filter({ r -> r.ref?.blast == null && r.ref?.fasta == null && r.ref?.sequence != null })
				.map({ r -> record(meta: r.meta, sequence: r.ref.sequence) })
		).map({ r -> record(meta: r.meta, ref_fasta: r.fasta) }))

	def genome_db_ch = MAKEBLASTDB(
		input_ch
			.join(ref_fasta_ch, by: 'meta')
			.map({ r -> record(fasta: record(meta: r.meta, path: r.genome.fasta), args: null) })
	)
	def blast_ch = BLASTN_REFERENCE(
		ref_fasta_ch
			.join(genome_db_ch.map({ d -> record(meta: d.meta, db: d) }), by: 'meta')
			.map({ r -> record(query: record(meta: r.meta, path: r.ref_fasta), subject_db: r.db, args: blastn_args) })
	).map({ r -> record(meta: r.meta_query, blast6: r.path) })
	def hits_ch = input_ch
		.filter({ r -> r.ref?.blast != null })
		.map({ r -> record(meta: r.meta, blast6: r.ref.blast) })
		.mix(blast_ch)

	// One output per genome: the rotated genome, or the input genome when it has no reference or no hit
	def genome_hits_ch = input_ch
		.join(hits_ch.map({ r -> record(meta: r.meta, hits: r.blast6) }), by: 'meta', remainder: true)
	def rotated_ch = ROTATE_FASTA(
		genome_hits_ch
			.filter({ r -> r.hits != null && !(r.hits as Path).isEmpty() })
			.map({ r -> record(meta: r.meta, fasta: r.genome.fasta, blast: r.hits) }),
		file("${moduleDir}/assets")
	)
		.map({ r -> record(meta: r.meta, fasta: r.fasta, gff: null, info: null) })    // the rotation invalidates the coordinates
		.mix(genome_hits_ch
			.filter({ r -> r.hits == null || (r.hits as Path).isEmpty() })
			.map({ r -> r.genome }))

	emit:
	rotated: Channel<Genome>       = rotated_ch
	hits: Channel<ReferenceHits>  = hits_ch
}
