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


// Annotate the FastANI results with the database information, in a folder with the database description and taxonomy
process FASTANI_TO_ORGFINDER {
	container 'registry.gitlab.unige.ch/amr-genomics/rscript:v3'
	label 'rscript'
	memory 8.GB
	cpus 1
	time 30.m

	input:
	record(meta: Map, fastani: Path, db_dir: Path)
	script_file: Path

	stage:
	stageAs(fastani, 'input.fastani')
	stageAs(db_dir, 'db')
	stageAs(script_file, 'fastani2orgfinder.R')

	output:
	record(
		meta: meta,
		dir: file('orgfinder'),            // ani.tsv, and db.tsv, tax.tsv, tax.rds of the database
		tsv: file('orgfinder/ani.tsv'),
		org_name: env('ORG_NAME'),    // best hit: organism name, assembly accession and ANI (empty when no hit)
		org_acc: env('ORG_ACC'),
		ani: env('ORG_ANI')
	)

	script:
	"""
	mkdir orgfinder
	Rscript fastani2orgfinder.R
	cp -L db/db.tsv db/tax.tsv db/tax.rds orgfinder/
	source org.env
	"""

	stub:
	"""
	mkdir orgfinder
	touch orgfinder/ani.tsv orgfinder/db.tsv orgfinder/tax.tsv orgfinder/tax.rds
	export ORG_NAME="xxx" ORG_ACC="xxx" ORG_ANI="90"
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
			.map({ r -> record(meta: r.meta, query: [r.fasta]) })
			.combine(refs: db.map({ d -> d.fasta }))
			.map({ r -> r + record(args: fastani_args) })
	)

	// Annotate ANI results with the database information (not available when references are given as FASTA)
	def orgfinder_ch = FASTANI_TO_ORGFINDER(
		fastani_ch
			.combine(db_dir: db.map({ d -> d.dir }))
			.filter({ r -> r.db_dir != null })
			.map({ r -> record(meta: r.meta, fastani: r.tsv, db_dir: r.db_dir) }),
		file("${moduleDir}/assets/fastani2orgfinder.R")
	)

	emit:
	fastani: Channel<Record>   = fastani_ch      // records with fields meta, tsv
	orgfinder: Channel<Record> = orgfinder_ch    // records with fields meta, dir, tsv, org_name, org_acc, ani
}
