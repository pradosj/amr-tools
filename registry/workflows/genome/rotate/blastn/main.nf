nextflow.enable.types = true

include { MAKEBLASTDB  } from './modules/blast/makeblastdb'
include { BLASTN       } from './modules/blast/blastn'
include { FASTA_ROTATE } from './modules/fasta/rotate'
include { MetaPath     } from './subworkflows/input/main.nf'


workflow GENOME_ROTATE_BLASTN {
	take:
	inputs: Channel<MetaPath>
	ref: Value<Path>    // FASTA of the reference dataset

	main:
	def fasta_ch = inputs.map({ r -> record(meta: r.meta, fasta: r) })
	def db = MAKEBLASTDB(ref.map({ p -> record(meta: null, path: p) }))    // BLAST database of the reference
	def hits_ch = BLASTN(
		inputs
			.combine(db: db)
			.map({ r -> record(meta: r.meta, path: r.path, ref: record(meta: null, path: r.db.path), args: '-perc_identity 70 -evalue 1e-3 -subject_besthit') })
	)

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
