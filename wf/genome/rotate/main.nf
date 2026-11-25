nextflow.enable.types = true

include { GENOME_ROTATE                        } from './workflows/rotate/main.nf'
include { INPUT_FASTA; MetaPath                } from './subworkflows/input/main.nf'
include { GZIP_UNZ                             } from './modules/io/gzip/unz'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'


params {
	fasta: String
	ref: Path
}


workflow {
	main:
	// Validate parameters and print summary of supplied ones
	validateParameters()
	log.info(paramsSummaryLog(workflow))

	// the reference dataset (decompressed when gzipped)
	def ref: Value<Path> = channel.value(params.ref)
	if (params.ref.name.endsWith('.gz')) {
		ref = GZIP_UNZ(channel.value(record(meta: record(name: 'reference'), path: params.ref))).map({ r -> r.file })
	}
	def genomes: Channel<MetaPath> = INPUT_FASTA(params.fasta)
	def rotated = GENOME_ROTATE(genomes.combine(ref: ref).map({ r -> record(fasta: record(meta: r.meta, path: r.path), origins: record(meta: record(id: 'reference'), path: r.ref)) }))

	publish:
	fasta = rotated.fasta
	rotation = rotated.rotation
}


output {
	fasta: Channel<MetaPath> {
		path { r -> r.path >> "genome-rotate/${r.meta.id}.fasta" }
	}
	rotation: Channel<MetaPath> {
		path { r -> r.path >> "genome-rotate/rotation/${r.meta.id}.tsv" }
	}
}
