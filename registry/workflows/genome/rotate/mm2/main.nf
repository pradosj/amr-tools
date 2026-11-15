nextflow.enable.types = true

include { MINIMAP2_ALIGN_PAF } from './modules/minimap2/align_paf'
include { FASTA_ROTATE       } from './modules/fasta/rotate'
include { MetaPath           } from './subworkflows/input/main.nf'

// A FASTA file to rotate on a reference dataset
record RotateInput {
	fasta: MetaPath    // the sequences to rotate (its meta, unique to each input, is the meta of the outputs)
	ref: MetaPath      // FASTA of the reference dataset
}


// Convert minimap2 PAF hits to the BLAST tabular format of FASTA_ROTATE (-outfmt "6 std qlen slen"): the number of
// matching bases as bit score, no e-value and gap opening (0), the mismatches being the aligned bases that do not match
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


// Rotate each sequence of FASTA files to start at its best hit (number of matching bases) on a reference dataset, on
// the strand of the hit: minimap2 of the FASTA (query) on the reference FASTA (default: -x asm5 -c, the base-level
// alignment giving the exact ends of the hits), then
// FASTA_ROTATE
workflow GENOME_ROTATE_MM2 {
	take:
	inputs: Channel<RotateInput>
	minimap2_args: String?

	main:
	def fasta_ch = inputs.map({ r -> record(meta: r.fasta.meta, fasta: r.fasta) })
	def paf_ch = MINIMAP2_ALIGN_PAF(inputs.map({ r -> record(meta: r.fasta.meta, path: r.fasta.path, ref: r.ref, args: minimap2_args ?: '-x asm5 -c') }))
	def hits_ch = PAF_TO_BLAST6(paf_ch)

	// One output per FASTA: the rotated FASTA, or the input FASTA when it has no hit
	def fasta_hits_ch = fasta_ch.join(hits_ch.map({ h -> record(meta: h.meta, hits: h) }), by: 'meta')
	def rotated_ch = FASTA_ROTATE(
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
