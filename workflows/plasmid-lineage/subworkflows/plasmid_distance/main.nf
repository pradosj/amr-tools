nextflow.enable.types = true

include { BLASTN as BLASTN_SELF } from './modules/blast/blastn'


// A set of sequences (e.g. plasmids): one FASTA file, one sequence per member of the set
record SequenceSet {
	meta: Record
	fasta: Path
}


// SNP + structural variant distance of each pair of sequences (plasmid_distance() of assets/lib_plasmid_distance.R),
// from their BLASTN against themselves
process SELF_BLASTN_DISTANCE_GRAPH {
	container 'registry.gitlab.unige.ch/amr-genomics/rscript:v3'
	label 'rscript'
	memory 4.GB
	cpus 1
	time 30.m

	input:
	record(meta: Record, blast: Path)
	assets: Path    // the assets folder of the subworkflow (lib_io.R, lib_plasmid_distance.R)

	stage:
	stageAs(blast, 'self.blastn')
	stageAs(assets, 'assets')

	output:
	record(
		meta: meta,
		graph: file('blast_distance_graph.rds'),     // distance matrix: column seq, then one column per sequence
	)

	script:
	"""
	#!/usr/bin/env Rscript
	source("assets/lib_io.R")
	source("assets/lib_plasmid_distance.R")
	g <- blast_distance_graph("self.blastn")
	saveRDS("blast_distance_graph.rds")
	"""

	stub:
	"""
	touch blast_distance_graph.rds
	"""
}


// SNP + structural variant distance matrix of the sequences of each set (e.g. the plasmids of a set, in one FASTA
// file): BLASTN of the sequences against themselves, then for each pair of sequences the mismatches, the gap openings
// and the structural variants (unaligned regions, breakpoints) of a chain of their alignments (see plasmid_distance()).
workflow PLASMID_DISTANCE {
	take:
	sets: Channel<SequenceSet>    // one FASTA file per set, one sequence per member of the set
	blastn_args: String?          // alignment of the sequences

	main:
	def blast_ch = BLASTN_SELF(sets.map({ r -> record(
		query: r,
		args: blastn_args ?: '-perc_identity 70 -evalue 1e-3 -subject_besthit'
	) }))
	def distance_ch = SELF_BLASTN_DISTANCE_GRAPH(
		blast_ch.map({ r -> record(meta: r.meta_query, blast: r.tsv) }),
		file("${moduleDir}/assets")
	)

	emit:
	distance: Channel<Record> = distance_ch    // records with fields meta (of the set), tsv (distance matrix), details
}
