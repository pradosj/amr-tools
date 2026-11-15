nextflow.enable.types = true

include { GENOME_ROTATE                        } from './workflows/rotate/main.nf'
include { INPUT_FASTA; MetaPath                } from './subworkflows/input/main.nf'
include { MAKEBLASTDB                          } from './modules/blast/makeblastdb'
include { GZIP_UNZ                             } from './modules/gzip/unz'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'


params {
	fasta: String
	ref: Path
	blastn_args: String = '-perc_identity 70 -evalue 1e-3 -subject_besthit'
}


workflow {
	main:
	// Validate parameters and print summary of supplied ones
	validateParameters()
	log.info(paramsSummaryLog(workflow))

	// BLAST database of the reference dataset (decompressed when gzipped)
	def ref: Value<Path> = channel.value(params.ref)
	if (params.ref.name.endsWith('.gz')) {
		ref = GZIP_UNZ(channel.value(record(meta: record(name: 'reference'), path: params.ref))).map({ r -> r.file })
	}
	def ref_db = MAKEBLASTDB(ref.map({ p -> record(meta: record(id: 'reference'), path: p) }))
	def genomes: Channel<MetaPath> = INPUT_FASTA(params.fasta)
	def rotated = GENOME_ROTATE(genomes.combine(ref: ref_db).map({ r -> record(fasta: record(meta: r.meta, path: r.path), ref: r.ref) }), params.blastn_args)

	publish:
	rotated_fasta = rotated.rotated_fasta
	hits = rotated.hits
}


output {
	rotated_fasta: Channel<MetaPath> {
		path { r -> r.path >> "genome-rotate/${r.meta.id}.fasta" }
	}
	hits: Channel<Record> {
		path { r -> r.path >> "genome-rotate/hits/${r.meta.id}.blastn" }
	}
}
