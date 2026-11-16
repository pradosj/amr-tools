nextflow.enable.types = true

include { MINIMAP2_ALIGN_BAM } from './modules/minimap2/align_bam'
include { FASTA_ROTATE_BAM   } from './modules/fasta/rotate'
include { MetaPath           } from './subworkflows/input/main.nf'



workflow GENOME_ROTATE_MM2 {
	take:
	inputs: Channel<MetaPath>
	ref: Value<Path>    // FASTA of the reference dataset

	main:
	def bam_ch = MINIMAP2_ALIGN_BAM(inputs.combine(ref: ref).map({ r -> record(meta: r.meta, path: r.path, ref: record(meta: null, path: r.ref), args: '-x asm5') })).bam

	// One output per FASTA (the sequences without alignment are not rotated)
	def rotated_ch = FASTA_ROTATE_BAM(
		inputs
			.map({ r -> record(meta: r.meta, fasta: r) })
			.join(bam_ch.map({ b -> record(meta: b.meta, bam: b) }), by: 'meta')
			.map({ r -> record(fasta: r.fasta, bam: r.bam) })
	)

	emit:
	rotated_fasta: Channel<MetaPath> = rotated_ch    // the rotated FASTA files (same meta)
	hits = bam_ch    // the minimap2 alignments of each FASTA (meta, meta_ref, path: BAM, the FASTA as query)
}
