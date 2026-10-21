
include { ORGFINDER_DB_BUILD } from './subworkflows/orgfinder_db_build/main.nf'
include { FASTANI } from './modules/fastani/main.nf'
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'


process FASTANI_TO_ORGFINDER {
	  label 'rscript'
	  container 'registry.gitlab.unige.ch/amr-genomics/rscript:v3'
    memory '8 GB'
    cpus 1
    time '30 min'
    input:
        tuple(val(meta), path('input.fastani'))
        each path('fastani2orgfinder.R')
        each path('db')
    output:
        tuple(val(meta), path('output.orgfinder'))
    script:
				"""
				Rscript fastani2orgfinder.R
				"""
}


workflow {
	main:
		// Validate parameters and print summary of supplied ones
		validateParameters()
		log.info(paramsSummaryLog(workflow))

		// Prepare reference database
		def ref_ch = Channel.empty()
		def db_dir = Channel.empty()
		if (params.ref.db != null) {
			ref_ch = Channel.fromPath("${params.ref.db}/fna/*.fna").collect().map({[it]})
		} else if (params.ref.fasta != null) {
			ref_ch = Channel.fromPath(params.ref.fasta).collect().map({[it]})
		} else {
			ORGFINDER_DB_BUILD()
			ref_ch = ORGFINDER_DB_BUILD.out.fasta.collect().map({[it]})
			db_dir = ORGFINDER_DB_BUILD.out.dir
		}

		// Prepare query genomes
		def query_ch = Channel.empty()
		if (params.query!=null) {
			query_ch = Channel.fromPath(params.query).map({
					def id = it.name.replaceAll(/\.(fasta|fna|fa)$/,'')
					[[sample_id:id],it]
			})
		}
		
		// Run FASTANI
		query_ch
				.combine(ref_ch)
				.map({m, q, r -> tuple(m, q, r)})
				| FASTANI
				
		FASTANI_TO_ORGFINDER(
			FASTANI.out,
			file("${moduleDir}/assets/fastani2orgfinder.R"),
			db_dir
		)

	publish:
		orgfinder_tsv = FASTANI_TO_ORGFINDER.out
		db = db_dir
}

output {
	db {
		path {x -> x >> "db/orgfinder"}
	}
	orgfinder_tsv {
		path { m,x -> x >> "${m.sample_id}.orgfinder"}
	}
}




