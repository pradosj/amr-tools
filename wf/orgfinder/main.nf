nextflow.enable.types = true

include { ORGFINDER; OrgfinderDb               } from './workflows/orgfinder/main.nf'
include { ORG_DB_BUILD                   } from './subworkflows/org_db/main.nf'
include { INPUT_FASTA as INPUT_REF_FASTA; INPUT_FASTA as INPUT_QUERY_FASTA } from './subworkflows/input/main.nf'
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


workflow {
	main:
	// Validate parameters and print summary of supplied ones
	validateParameters()
	log.info(paramsSummaryLog(workflow))

	// Reference database: given (--ref.db), FASTA files (--ref.fasta), or built from NCBI genomes
	def db_dir_ch = channel.empty()
	def db: Value<OrgfinderDb>
	if (params.ref?.db != null) {
		db = channel.value(record(dir: file(params.ref.db), fasta: files("${params.ref.db}/fna/*.fna").toSorted(), tsv: file("${params.ref.db}/db.tsv")))
	} else if (params.ref?.fasta != null) {
		db = INPUT_REF_FASTA(params.ref.fasta).map({ r -> r.path }).collect().map({ fs -> record(fasta: fs.toList().toSorted({ f -> f.name })) })
	} else {
		db = ORG_DB_BUILD()
		db_dir_ch = channel.empty().mix(db.map({ d -> d.dir }))
	}

	// Genomes to identify
	def query_ch = INPUT_QUERY_FASTA(params.query).map({ r -> record(meta: record(sample_id: r.meta.id), fasta: r.path) })    // gzipped files are decompressed

	def orgfinder = ORGFINDER(db, query_ch, params.fastani_args)

	publish:
	fastani_tsv = orgfinder.fastani
	orgfinder_tsv = orgfinder.orgfinder
	db = db_dir_ch
}

output {
	db: Channel<Path> {
		path { x -> x >> "db/org_db" }
	}
	orgfinder_tsv: Channel<Record> {
		path { r -> r.tsv >> "${r.meta.sample_id}.orgfinder" }
	}
	fastani_tsv: Channel<Record> {
		path { r -> r.tsv >> "${r.meta.sample_id}.fastani" }
	}
}
