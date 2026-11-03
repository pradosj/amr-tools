nextflow.enable.types = true

include { FASTANI } from './modules/fastani/main.nf'


// Genome to identify
record Genome {
	meta: Map
	fasta: Path
	gff: Path?
	info: Path?
}

// Reference database: FASTA of the reference genomes, with the database folder and its description table (db.tsv)
// when the database was built by ORG_DB_BUILD (dir and tsv are null for references given as FASTA files)
record OrgfinderDb {
	dir: Path?
	fasta: List<Path>
	tsv: Path?
}


// Annotate the FastANI results with the database information
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


// Identify the organism of each genome by its ANI with the reference genomes of the database
workflow ORGFINDER {
	take:
	db: Value<OrgfinderDb>
	genomes: Channel<Genome>
	fastani_args: String?

	main:
	// Compute ANI of the genomes against all reference genomes
	def fastani_ch = FASTANI(
		genomes
			.map({ r -> record(meta: r.meta, query: r.fasta) })
			.combine(refs: db.map({ d -> d.fasta }))
			.map({ r -> r + record(args: fastani_args) })
	)

	// Annotate ANI results with the database information (not available when references are given as FASTA)
	def orgfinder_ch = FASTANI_TO_ORGFINDER(
		fastani_ch
			.combine(db_tsv: db.map({ d -> d.tsv }))
			.filter({ r -> r.db_tsv != null })
			.map({ r -> record(meta: r.meta, fastani: r.tsv, db_tsv: r.db_tsv) }),
		file("${moduleDir}/assets/fastani2orgfinder.R")
	)

	emit:
	fastani: Channel<Record>   = fastani_ch      // records with fields meta, tsv
	orgfinder: Channel<Record> = orgfinder_ch    // records with fields meta, tsv
}
