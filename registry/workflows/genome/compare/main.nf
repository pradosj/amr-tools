nextflow.enable.types = true

include { MINIMAP2_ALIGN_BAM                        } from './modules/minimap2/align_bam'
include { MINIMAP2_ALIGN_BAM as ROTATION_ALIGNMENT  } from './modules/minimap2/align_bam'
include { BCFTOOLS_ASM5_MPILEUP                     } from './modules/bcftools/asm5_mpileup'
include { QUARTO_RENDER                             } from './modules/quarto/render'


// A query assembly to compare with the reference assembly
record Comparison {
	meta: Record
	ref: Path
	query: Path
}


process ROTATE_FASTA {
	container "registry.gitlab.unige.ch/amr-genomics/rscript:v3"
	label 'rscript'
	memory 8.GB
	cpus 1
	time 30.m

	input:
	record(meta: Record, query: Path, bam: Path)
	assets: Path

	stage:
	stageAs(query, 'query.fasta')
	stageAs(bam, 'aln.bam')
	stageAs(assets, 'assets')

	output:
	record(meta: meta, fasta: file('rotated_query.fasta'))

	script:
	"""
	#!/usr/bin/env Rscript
	source("assets/lib_bam.R")
	rotate_queries("query.fasta","aln.bam") |>
		writeXStringSet("rotated_query.fasta")
	"""

	stub:
	"""
	touch rotated_query.fasta
	"""
}

// Rotate the contigs of the queries to start at the same position as in the reference
workflow ROTATE_QUERIES {
	take:
	queries: Channel<Comparison>
	minimap2_args: String?

	main:
	def aln_ch = ROTATION_ALIGNMENT(queries.map({ r -> record(meta: r.meta, path: r.query, ref: record(meta: null, path: r.ref), args: "-x asm5 ${minimap2_args ?: ''}") })).bam
	def rotated_ch = ROTATE_FASTA(
		queries
			.join(aln_ch, by: 'meta')
			.map({ r -> record(meta: r.meta, query: r.query, bam: r.path) }),
		file("${moduleDir}/assets")
	)
	def out_ch = queries
		.join(rotated_ch, by: 'meta')
		.map({ r -> record(meta: r.meta, ref: r.ref, query: r.fasta) })

	emit:
	rotated_queries: Channel<Comparison> = out_ch
}


// Compare query assemblies with their reference assembly: rotation of the queries on the reference (unless
// no_query_rotation), alignment (minimap2 -x asm5), mutations (bcftools mpileup) and HTML report
workflow GENOME_COMPARE {
	take:
	queries: Channel<Comparison>
	ref_gff: Path?                // annotation of the reference, for the consequences of the mutations
	no_query_rotation: Boolean
	minimap2_args: String?
	mpileup_args: String?

	main:
	def query_ch: Channel<Comparison> = queries
	def rotated_query_ch: Channel<Comparison> = channel.empty()
	if (!no_query_rotation) {
		query_ch = ROTATE_QUERIES(queries, minimap2_args)
		rotated_query_ch = query_ch
	}

	// Map (rotated) queries to their reference
	def aln_ch = MINIMAP2_ALIGN_BAM(
		query_ch.map({ r -> record(meta: r.meta, path: r.query, ref: record(meta: null, path: r.ref), args: "-x asm5 ${minimap2_args ?: ''}") })
	).outputs

	// Call mutations
	def mut_ch = BCFTOOLS_ASM5_MPILEUP(
		aln_ch
			.join(query_ch.map({ r -> record(meta: r.meta, ref: r.ref) }), by: 'meta')
			.map({ r -> record(meta: r.meta, ref: r.ref, bam: r.bam, gff: ref_gff, args: mpileup_args) })
	)

	// HTML report
	def html_ch = QUARTO_RENDER(
		aln_ch
			.join(mut_ch, by: 'meta')
			.map({ r -> record(
				meta: r.meta,
				template: file("${moduleDir}/assets/bam_report.qmd"),
				extra: file("${moduleDir}/assets"),
				files: [r.bam, r.txt],
				args: "-P 'bam_file:${(r.bam as Path).name}' -P 'txt_file:${(r.txt as Path).name}'"
			) })
	)

	emit:
	rotated_queries = rotated_query_ch    // records with fields meta, ref, query (rotated query)
	alignment = aln_ch                    // records with fields meta, meta_ref, bam, bai
	mutations = mut_ch                    // records with fields meta, vcf, csi, txt
	html = html_ch                        // records with fields meta, html
}
