nextflow.enable.types = true

include { MINIMAP2_ALIGN_BAM } from './modules/aln/minimap2/align_bam'
include { FASTA_ROTATE_BAM   } from './modules/fa/rotate'
include { MetaPath           } from './subworkflows/input/main.nf'

// A FASTA file to rotate on the origins of its sequences
record RotateInput {
	fasta: MetaPath      // the sequences to rotate (its meta, unique to each input, is the meta of the outputs)
	origins: MetaPath    // FASTA of the reference sequences to start at
}


workflow GENOME_ROTATE {
	take:
	inputs: Channel<RotateInput>

	main:
	def bam_ch = MINIMAP2_ALIGN_BAM(inputs.map({ r -> record(meta: r.fasta.meta, path: r.fasta.path, ref: r.origins, args: '-x asm5') })).bam

	// One output per FASTA (the sequences without alignment are not rotated)
	def rotated_ch = FASTA_ROTATE_BAM(
		inputs
			.map({ r -> record(meta: r.fasta.meta, fasta: r.fasta) })
			.join(bam_ch.map({ b -> record(meta: b.meta, bam: b) }), by: 'meta')
			.map({ r -> record(fasta: r.fasta, bam: r.bam) })
	)

	emit:
	fasta: Channel<MetaPath> = rotated_ch.fasta    // the rotated FASTA files (same meta)
	rotation: Channel<MetaPath> = rotated_ch.tsv    // the rotation of each FASTA (TSV: chr, pos, strand of the new start of each sequence)
}
