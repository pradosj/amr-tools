nextflow.enable.types = true

include { ORGFINDER_DB_BUILD } from './subworkflows/orgfinder_db_build/main.nf'
include { FASTANI            } from './modules/fastani/main.nf'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

// Reference genomes: either FASTA files (glob) or a directory with a database previously built by the pipeline
record RefParams {
	fasta: String?
	db: String?
}

params {
	query: String? = null
	ref: RefParams = record(fasta: null, db: null)
	fastani_args: String? = null
}


process FASTANI_TO_ORGFINDER {
	container 'registry.gitlab.unige.ch/amr-genomics/rscript:v3'
	label 'rscript'
	memory 8.GB
	cpus 1
	time 30.m

	input:
	record(meta: Map, fastani: Path, db_tsv: Path)
	script_file: Path

	stage:
	stageAs(fastani, 'input.fastani')
	stageAs(db_tsv, 'db.tsv')
	stageAs(script_file, 'fastani2orgfinder.R')

	output:
	record(meta: meta, tsv: file('output.orgfinder'))

	script:
	"""
	Rscript fastani2orgfinder.R
	"""

	stub:
	"""
	touch output.orgfinder
	"""
}


workflow {
	main:
	// Validate parameters and print summary of supplied ones
	validateParameters()
	log.info(paramsSummaryLog(workflow))

	// Prepare reference database
	def db_dir_ch = channel.empty()
	def refs: Value<List<Path>>
	def db_tsv: Value<Path?>
	if (params.ref?.db != null) {
		refs = channel.value(files("${params.ref.db}/fna/*.fna").toSorted())
		db_tsv = channel.value(file("${params.ref.db}/db.tsv"))
	} else if (params.ref?.fasta != null) {
		refs = channel.value(files(params.ref.fasta).toSorted())
		db_tsv = channel.value(null)
	} else {
		def db = ORGFINDER_DB_BUILD()
		refs = db.map({ d -> d.fasta })
		db_tsv = db.map({ d -> d.tsv })
		db_dir_ch = channel.empty().mix(db.map({ d -> d.dir }))
	}

	// Prepare query genomes
	def query_ch = channel.empty()
	if (params.query != null) {
		query_ch = channel.fromPath(params.query).map({ f ->
			def id = f.name.replaceAll(/\.(fasta|fna|fa)$/, '')
			record(meta: [sample_id: id], query: f)
		})
	}

	// Compute ANI of the queries against all reference genomes
	def fastani_ch = FASTANI(
		query_ch
			.combine(refs: refs)
			.map({ r -> r + record(args: params.fastani_args) })
	)

	// Annotate ANI results with the database information (not available when references are given as FASTA)
	def orgfinder_ch = FASTANI_TO_ORGFINDER(
		fastani_ch
			.combine(db_tsv: db_tsv)
			.filter({ r -> r.db_tsv != null })
			.map({ r -> record(meta: r.meta, fastani: r.tsv, db_tsv: r.db_tsv) }),
		file("${projectDir}/assets/fastani2orgfinder.R")
	)

	publish:
	fastani_tsv = fastani_ch
	orgfinder_tsv = orgfinder_ch
	db = db_dir_ch
}

output {
	db: Channel<Path> {
		path { x -> x >> "db/orgfinder" }
	}
	orgfinder_tsv: Channel<Record> {
		path { r -> r.tsv >> "${r.meta.sample_id}.orgfinder" }
	}
	fastani_tsv: Channel<Record> {
		path { r -> r.tsv >> "${r.meta.sample_id}.fastani" }
	}
}
