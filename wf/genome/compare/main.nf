nextflow.enable.types = true

include { GENOME_COMPARE                       } from './workflows/compare/main.nf'
include { INPUT_FASTA; MetaPath                } from './subworkflows/input/main.nf'
include { GZIP_UNZ                             } from './modules/io/gzip/unz'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

params {
	ref: Path
	ref_gff: Path? = null
	query: String
	no_query_rotation: Boolean = false
	minimap2_args: String? = '-z1000000 -p 0.1 -Y --cs'
	mpileup_args: String? = '--no-BAQ --indel-size 1000 --ff=UNMAP --ambig-reads=drop --gap-frac=0 --min-ireads=0'
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
	def query_ch: Channel<MetaPath> = INPUT_FASTA(params.query)
		.map({ r -> record(meta: record(id: r.meta.id), path: r.path) })

	def cmp = GENOME_COMPARE(query_ch, ref, params.ref_gff, params.no_query_rotation, params.minimap2_args, params.mpileup_args)

	publish:
	rotated_query_fa = cmp.rotated_queries
	alignment = cmp.alignment
	mutations = cmp.mutations
	html = cmp.html
}

output {
	rotated_query_fa: Channel<MetaPath> {
		path { r -> r.path >> "${r.meta.id}.rotated.fasta" }
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
