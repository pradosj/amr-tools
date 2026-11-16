nextflow.enable.types = true

include { GENOME_ROTATE         } from './workflows/rotate/main.nf'
include { MINIMAP2_ALIGN_BAM    } from './modules/minimap2/align_bam'
include { BCFTOOLS_ASM5_MPILEUP } from './modules/bcftools/asm5_mpileup'
include { QUARTO_RENDER         } from './modules/quarto/render'
include { MetaPath              } from './subworkflows/input/main.nf'


// Compare query assemblies with a reference assembly: rotation of the queries on the reference with GENOME_ROTATE
// (unless no_query_rotation), alignment (minimap2 -x asm5), mutations (bcftools mpileup) and HTML report
workflow GENOME_COMPARE {
	take:
	queries: Channel<MetaPath>    // the FASTA of each query assembly (meta unique to each one)
	ref: Value<Path>              // FASTA of the reference assembly
	ref_gff: Path?                // annotation of the reference, for the consequences of the mutations
	no_query_rotation: Boolean
	minimap2_args: String?
	mpileup_args: String?

	main:
	def query_ch: Channel<MetaPath> = queries
	def rotated_query_ch: Channel<MetaPath> = channel.empty()
	if (!no_query_rotation) {
		query_ch = GENOME_ROTATE(queries.combine(ref: ref).map({ r -> record(fasta: record(meta: r.meta, path: r.path), origins: record(meta: null, path: r.ref)) })).fasta
		rotated_query_ch = query_ch
	}

	// Map (rotated) queries to the reference
	def aln_ch = MINIMAP2_ALIGN_BAM(
		query_ch.combine(ref: ref).map({ r -> record(meta: r.meta, path: r.path, ref: record(meta: null, path: r.ref), args: "-x asm5 ${minimap2_args ?: ''}") })
	).outputs

	// Call mutations
	def mut_ch = BCFTOOLS_ASM5_MPILEUP(
		aln_ch.combine(ref: ref).map({ r -> record(meta: r.meta, ref: r.ref, bam: r.bam, gff: ref_gff, args: mpileup_args) })
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
	rotated_queries = rotated_query_ch    // the rotated queries (meta, path)
	alignment = aln_ch                    // records with fields meta, meta_ref, bam, bai
	mutations = mut_ch                    // records with fields meta, vcf, csi, txt
	html = html_ch                        // records with fields meta, html
}
