nextflow.enable.types = true

include { MINIMAP2_ALIGN_PAF } from './modules/minimap2/align_paf'
include { FASTA_ROTATE_BLAST } from './modules/fasta/rotate'
include { MetaPath           } from './subworkflows/input/main.nf'



process PAF_TO_BLAST6 {
	container 'docker.io/ubuntu:24.04'
	label 'local'
	memory 1.GB
	cpus 1
	time 30.m

	input:
	record(meta: Record?, meta_ref: Record?, path: Path)

	stage:
	stageAs(path, 'hits.paf')

	output:
	record(meta: meta, meta_ref: meta_ref, path: file('hits.blast6'))

	script:
	"""
	awk 'BEGIN { FS = OFS = "\\t" } {
		plus = (\$5 == "+")
		print \$1, \$6, sprintf("%.3f", 100 * \$10 / \$11), \$11, \$11 - \$10, 0, \$3 + 1, \$4, (plus ? \$8 + 1 : \$9), (plus ? \$9 : \$8 + 1), 0, \$10, \$2, \$7
	}' hits.paf > hits.blast6
	"""

	stub:
	"""
	touch hits.blast6
	"""
}


workflow GENOME_ROTATE_MM2 {
	take:
	inputs: Channel<MetaPath>
	ref: Value<Path>    // FASTA of the reference dataset

	main:
	def fasta_ch = inputs.map({ r -> record(meta: r.meta, fasta: r) })
	def paf_ch = MINIMAP2_ALIGN_PAF(inputs.combine(ref: ref).map({ r -> record(meta: r.meta, path: r.path, ref: record(meta: null, path: r.ref), args: '-x asm5 -c') }))
	def hits_ch = PAF_TO_BLAST6(paf_ch)

	// One output per FASTA: the rotated FASTA, or the input FASTA when it has no hit
	def fasta_hits_ch = fasta_ch.join(hits_ch.map({ h -> record(meta: h.meta, hits: h) }), by: 'meta')
	def rotated_ch = FASTA_ROTATE_BLAST(
		fasta_hits_ch
			.filter({ r -> !(r.hits.path as Path).isEmpty() })
			.map({ r -> record(fasta: r.fasta, blast6: r.hits) })
	)
		.mix(fasta_hits_ch
			.filter({ r -> (r.hits.path as Path).isEmpty() })
			.map({ r -> r.fasta }))

	emit:
	rotated_fasta: Channel<MetaPath> = rotated_ch    // the rotated FASTA files (same meta)
	hits = paf_ch    // the minimap2 hits of each FASTA (meta, meta_ref, path: PAF, the FASTA as query)
}
