
include { ORGFINDER_DB } from './subworkflows/orgfinder_db/main.nf'
include { FASTANI } from './modules/fastani/main.nf'
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

workflow {
	main:
		// Validate parameters and print summary of supplied ones
		validateParameters()
		log.info(paramsSummaryLog(workflow))


		def ref_ch = Channel.empty()
		if (params.ref_db != null) {
			ref_ch = Channel.fromPath("${params.ref_db}/fna/*.fna").collect().map({[it]})
		} else if (params.ref_fasta != null) {
			ref_ch = Channel.fromPath(params.ref_fasta).collect().map({[it]})
		} else {
			ORGFINDER_DB()
			ref_ch = ORGFINDER_DB.out.fasta.collect().map({[it]})
		}

		def query_ch = Channel.fromPath(params.query).map({
				def id = it.name.replaceAll(/\.(fasta|fna|fa)$/,'')
				[[sample_id:id],it]
		})

		query_ch
				.combine(ref_ch)
				.map({m, q, r -> tuple(m, q, r)})
				| FASTANI

	publish:
		fastani_tsv = FASTANI.out
		db = ORGFINDER_DB.out.dir
}

output {
	db {
		path {x -> x >> "db/orgfinder"}
	}
	fastani_tsv {
		path { m,x -> x >> "${m.sample_id}.fastani"}
	}
}




