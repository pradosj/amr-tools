nextflow.enable.types = true

include { MINIMAP2_ALIGN_ASM5                       } from './modules/minimap2/align_asm5'
include { MINIMAP2_ALIGN_ASM5 as ROTATION_ALIGNMENT } from './modules/minimap2/align_asm5'
include { BCFTOOLS_ASM5_MPILEUP                     } from './modules/bcftools/asm5_mpileup'
include { QUARTO_RENDER                             } from './modules/quarto/render'
include { INPUT_FASTA                               } from './subworkflows/input/main.nf'
include { GZIP_UNZ                                  } from './modules/gzip/unz'
include { validateParameters; paramsSummaryLog      } from 'plugin/nf-schema'

params {
	ref: Path
	ref_gff: Path? = null
	query: String
	no_query_rotation: Boolean = false
	minimap2_args: String? = '-z1000000 -p 0.1 -Y --cs -a'
	mpileup_args: String? = '--no-BAQ --indel-size 1000 --ff=UNMAP --ambig-reads=drop --gap-frac=0 --min-ireads=0'
}

// A query assembly to compare with the reference assembly
record Comparison {
	meta: Record
	ref: Path
	query: Path
}


process ROTATE_FASTA {
	container "registry.gitlab.unige.ch/amr-genomics/rscript:v2"
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
	def aln_ch = ROTATION_ALIGNMENT(queries.map({ r -> record(meta: r.meta, path: r.query, ref: record(meta: null, path: r.ref), args: minimap2_args) })).bam
	def rotated_ch = ROTATE_FASTA(
		queries
			.join(aln_ch, by: 'meta')
			.map({ r -> record(meta: r.meta, query: r.query, bam: r.path) }),
		file("${projectDir}/assets")
	)
	def out_ch = queries
		.join(rotated_ch, by: 'meta')
		.map({ r -> record(meta: r.meta, ref: r.ref, query: r.fasta) })

	emit:
	rotated_queries: Channel<Comparison> = out_ch
}


workflow {
	main:
	// Validate parameters and print summary of supplied ones
	validateParameters()
	log.info(paramsSummaryLog(workflow))

	// Reference and queries: gzipped files are decompressed
	def ref: Value<Path> = channel.value(params.ref)
	if (params.ref.name.endsWith('.gz')) {
		ref = GZIP_UNZ(channel.value(record(meta: record(name: 'reference'), path: params.ref))).map({ r -> r.file })
	}
	def query_ch: Channel<Comparison> = INPUT_FASTA(params.query)
		.combine(ref: ref)
		.map({ r -> record(meta: record(id: r.meta.id), ref: r.ref, query: r.path) })

	def rotated_query_ch = channel.empty()
	if (!params.no_query_rotation) {
		query_ch = ROTATE_QUERIES(query_ch, params.minimap2_args)
		rotated_query_ch = query_ch
	}

	// Map (rotated) queries to REF
	def aln_ch = MINIMAP2_ALIGN_ASM5(
		query_ch.map({ r -> record(meta: r.meta, path: r.query, ref: record(meta: null, path: r.ref), args: params.minimap2_args) })
	).outputs

	// Call mutations
	def mut_ch = BCFTOOLS_ASM5_MPILEUP(
		aln_ch.combine(ref: ref).map({ r -> record(meta: r.meta, ref: r.ref, bam: r.bam, gff: params.ref_gff, args: params.mpileup_args) })
	)

	// HTML report
	def html_ch = QUARTO_RENDER(
		aln_ch
			.join(mut_ch, by: 'meta')
			.map({ r -> record(
				meta: r.meta,
				template: file("${projectDir}/assets/bam_report.qmd"),
				extra: file("${projectDir}/assets"),
				files: [r.bam, r.txt],
				args: "-P 'bam_file:${(r.bam as Path).name}' -P 'txt_file:${(r.txt as Path).name}'"
			) })
	)

	publish:
	rotated_query_fa = rotated_query_ch
	alignment = aln_ch
	mutations = mut_ch
	html = html_ch
}

output {
	rotated_query_fa: Channel<Comparison> {
		path { r -> r.query >> "${r.meta.id}.rotated.fasta" }
	}
	alignment: Channel<Record> {
		path { r ->
			r.bam >> "${r.meta.id}.bam"
			r.bai >> "${r.meta.id}.bai"
		}
	}
	mutations: Channel<Record> {
		path { r ->
			r.vcf >> "${r.meta.id}.vcf.gz"
			r.csi >> "${r.meta.id}.vcf.gz.csi"
			r.txt >> "${r.meta.id}.txt"
		}
	}
	html: Channel<Record> {
		path { r -> r.html >> "${r.meta.id}.html" }
	}
}
