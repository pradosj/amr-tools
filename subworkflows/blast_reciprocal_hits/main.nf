nextflow.enable.types = true

// The BLASTP of the sequences of a genome (meta1, query) against those of a genome (meta2, subject), in the format
// read by read_blast6() (-outfmt "6 std qlen slen"). The metas name the genomes by their sample_id.
record PairBlast {
	meta1: Record
	meta2: Record
	blastp: Path
}


// Reciprocal hits of a pair of genomes, from its two BLASTP (blast2 null for a genome with itself: blast1 is used
// for both directions, each pair of distinct sequences being kept once)
process RECIPROCAL_HITS {
	container 'registry.gitlab.unige.ch/amr-genomics/rscript:v3'
	label 'rscript'
	memory 4.GB
	cpus 1
	time 30.m

	input:
	record(meta1: Record, meta2: Record, blast1: Path, blast2: Path?, min_cover: Float?)
	libs: List<Path>    // lib_io.R, lib_blast.R

	stage:
	stageAs(blast1, 'blast1.tsv')
	stageAs(blast2, 'blast2.tsv')
	stageAs(libs, 'lib/*')

	output:
	record(meta1: meta1, meta2: meta2, tsv: file('reciprocal_hits.tsv'))

	script:
	"""
	#!/usr/bin/env Rscript
	source("lib/lib_io.R")
	source("lib/lib_blast.R")
	hits <- blast_reciprocal_hits("blast1.tsv", "${blast2 ? 'blast2.tsv' : 'blast1.tsv'}"${min_cover != null ? ", min_cover = ${min_cover}" : ''})
	${blast2 ? '' : 'hits <- filter(hits, query_id < subject_id)'}
	hits |>
		mutate(query_genome = '${meta1.sample_id}', subject_genome = '${meta2.sample_id}', .before = 1) |>
		write_tsv("reciprocal_hits.tsv")
	"""

	stub:
	"""
	touch reciprocal_hits.tsv
	"""
}


// Reciprocal BLASTP hits of genomes: pairs the BLASTP of each pair of genomes (A against B with B against A, and a
// genome against itself), runs blast_reciprocal_hits() (R/lib_blast.R) on each pair, and gives the hits of each pair
// in a TSV: query_genome, subject_genome (sample_id of the metas), query_id, subject_id (ids of the sequences in their
// genome), pct_ident, align_len, evalue and bit_score of the worst direction. Each pair of genomes is given once; for a
// genome with itself, each pair of distinct sequences once.
workflow BLAST_RECIPROCAL_HITS {
	take:
	blastp: Channel<PairBlast>    // the BLASTP of each pair of genomes, in both directions
	min_cover: Float?             // minimum coverage of the query protein by its alignment (null: 0.75)

	main:
	// Each pair of genomes once (meta1 <= meta2 by sample_id), with the BLASTP of its two directions
	def pairs_ch = blastp
		.map({ r -> record(pair: record(meta1: r.meta1, meta2: r.meta2), blast1: r.blastp) })
		.join(blastp.map({ r -> record(pair: record(meta1: r.meta2, meta2: r.meta1), blast2: r.blastp) }), by: 'pair')
		.filter({ r -> (r.pair.meta1.sample_id as String) <= (r.pair.meta2.sample_id as String) })
		.map({ r -> record(
			meta1: r.pair.meta1,
			meta2: r.pair.meta2,
			blast1: r.blast1,
			blast2: r.pair.meta1 == r.pair.meta2 ? null : r.blast2,    // a genome with itself: a single BLASTP
			min_cover: min_cover
		) })
	def hits_ch = RECIPROCAL_HITS(pairs_ch, [file("${moduleDir}/assets/lib_io.R"), file("${moduleDir}/assets/lib_blast.R")])

	emit:
	hits: Channel<Record> = hits_ch    // records with fields meta1, meta2, tsv (the reciprocal hits of the pair)
}
