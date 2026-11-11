nextflow.enable.types = true


record PairBlast {
	meta_query: Record
	meta_subject: Record
	blast6: Path
}


// Reciprocal hits of a pair of genomes, from its two BLAST (blast_rev null for a genome with itself: blast_fwd is then
// used for both directions; the same file given twice would be staged only once, under one of the two names)
process RECIPROCAL_HITS {
	container 'registry.gitlab.unige.ch/amr-genomics/rscript:v3'
	label 'rscript'
	memory 4.GB
	cpus 1
	time 30.m

	input:
	record(meta_query: Record, meta_subject: Record, blast_fwd: Path, blast_rev: Path?, min_cover: Float)
	assets: Path    // the assets folder of the workflow (lib_io.R, lib_blast.R)

	stage:
	stageAs(blast_fwd, 'blast_fwd.tsv')    // the BLAST outputs share the same file name
	stageAs(blast_rev, 'blast_rev.tsv')
	stageAs(assets, 'assets')

	output:
	record(meta_query: meta_query, meta_subject: meta_subject, blast6: file('reciprocal_hits.tsv'))

	script:
	"""
	#!/usr/bin/env Rscript
	source("assets/lib_io.R")
	source("assets/lib_blast.R")
	hits <- blast_reciprocal_hits("blast_fwd.tsv", "${blast_rev ? 'blast_rev.tsv' : 'blast_fwd.tsv'}", min_cover = ${min_cover})
	write_tsv(hits, "reciprocal_hits.tsv", col_names = FALSE)
	"""

	stub:
	"""
	touch reciprocal_hits.tsv
	"""
}


// Reciprocal BLAST hits of genomes: pairs the BLAST (BLASTP or BLASTN) of each pair of genomes (A against B with B
// against A, and a genome against itself), runs blast_reciprocal_hits() (R/lib_blast.R) on each pair, and gives the
// hits of each pair (scores of the worst direction) in the BLAST tabular format (-outfmt "6 std qlen slen", read by
// read_blast6()), the query in the genome of meta_query and the subject in that of meta_subject. Each ordered pair of
// genomes is given (A-B and B-A, the same hits mirrored), a genome with itself included (each sequence with itself and
// each pair of distinct sequences in both orders).
workflow BLAST_RECIPROCAL_HITS {
	take:
	blast: Channel<PairBlast>    // the BLAST of each pair of genomes, in both directions (e.g. the output of BLASTP)
	min_cover: Float             // minimum coverage of the query sequence by its alignment (e.g. 0.75)

	main:
	// Each ordered pair of genomes, with the BLAST of its two directions
	def pairs_ch = blast
		.map({ r -> record(pair: record(meta_query: r.meta_query, meta_subject: r.meta_subject), blast_fwd: r.blast6) })
		.join(
			blast.map({ r -> record(pair: record(meta_query: r.meta_subject, meta_subject: r.meta_query), blast_rev: r.blast6) }),
			by: 'pair'
		)
		.map({ r -> record(
			meta_query: r.pair.meta_query,
			meta_subject: r.pair.meta_subject,
			blast_fwd: r.blast_fwd,
			blast_rev: r.pair.meta_query == r.pair.meta_subject ? null : r.blast_rev,    // a genome with itself: a single BLAST
			min_cover: min_cover
		) })
	def hits_ch = RECIPROCAL_HITS(pairs_ch, file("${moduleDir}/assets"))

	emit:
	hits: Channel<PairBlast> = hits_ch    // the reciprocal hits of each pair (blast6 format, as the BLAST)
}

