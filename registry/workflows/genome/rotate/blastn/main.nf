nextflow.enable.types = true

include { BLASTN       } from './modules/blast/blastn'
include { FASTA_ROTATE } from './modules/fasta/rotate'
include { MetaPath     } from './subworkflows/input/main.nf'

// A FASTA file to rotate on a reference dataset
record RotateInput {
	fasta: MetaPath    // the sequences to rotate (its meta, unique to each input, is the meta of the outputs)
	ref: MetaPath      // BLAST database of the reference dataset (output of MAKEBLASTDB)
}


// Rotate each sequence of FASTA files to start at its best hit (bit score) on a reference dataset, on the strand of
// the hit: BLASTN of the FASTA (query) against the reference database (subject), then FASTA_ROTATE
workflow GENOME_ROTATE {
	take:
	inputs: Channel<RotateInput>
	blastn_args: String?

	main:
	def fasta_ch = inputs.map({ r -> record(meta: r.fasta.meta, fasta: r.fasta) })
	def hits_ch = BLASTN(inputs.map({ r -> record(meta: r.fasta.meta, path: r.fasta.path, ref: r.ref, args: blastn_args) }))

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
	hits = hits_ch    // the hits of each FASTA (meta, meta_ref, path: -outfmt "6 std qlen slen", the FASTA as query)
}
