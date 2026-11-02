nextflow.enable.types = true

include { NCBI_DATASET_DOWNLOAD_GENOME } from './modules/ncbi/dataset/main.nf'
include { UNTAR                        } from './modules/tar/untar.nf'

// Reference database built from NCBI genomes
record OrgfinderDb {
	dir: Path
	fasta: List<Path>
	tsv: Path
}

// Make the tsv file with all accession numbers
process STEP1_ACCESSIONS_TSV {
	container 'docker.io/staphb/ncbi-datasets:18.18.0'
	memory 8.GB
	cpus 1
	time 30.m

	input:
	datasets: Bag<Path>

	stage:
	stageAs(datasets, 'genomes/query*')

	output:
	record(genomes: file('genomes'), accessions: file('accessions.tsv'))

	script:
	"""
	cat genomes/query*/ncbi_dataset/data/assembly_data_report.jsonl \\
		| dataformat tsv genome --force --fields accession,organism-name,organism-tax-id,assmstats-total-sequence-len,assmstats-total-number-of-chromosomes \\
		| sort \\
		| uniq \\
		> accessions.tsv
	"""

	stub:
	"""
	touch accessions.tsv
	"""
}

process STEP2_MOLECULE_TYPES_TSV {
	container 'registry.gitlab.unige.ch/amr-genomics/rscript:v3'
	memory 2.GB
	cpus 1
	time 30.m

	input:
	record(genomes: Path, accessions: Path)

	stage:
	stageAs(genomes, 'ncbi_db/genomes')
	stageAs(accessions, 'ncbi_db/accessions.tsv')

	output:
	file('ncbi_db')

	script:
	"""
	jq -r '[.assemblyAccession, .chrName, .assignedMoleculeLocationType] | @tsv' ncbi_db/genomes/query*/ncbi_dataset/data/*/sequence_report.jsonl \\
		| sort \\
		| uniq \\
		> ncbi_db/molecule_types.tsv
	"""

	stub:
	"""
	touch ncbi_db/molecule_types.tsv
	"""
}

process STEP3_DB_BUILD {
	container 'registry.gitlab.unige.ch/amr-genomics/rscript:v3'
	memory 8.GB
	cpus 1
	time 30.m

	input:
	ncbi_db: Path
	taxdump: Path
	db_build_script: Path

	stage:
	stageAs(ncbi_db, 'ncbi_db')
	stageAs(taxdump, 'taxdump')
	stageAs(db_build_script, 'db_build.R')

	output:
	record(
		dir: file('orgfinder_db'),
		fasta: files('orgfinder_db/fna/*.fna').toSorted(),
		tsv: file('orgfinder_db/db.tsv')
	)

	script:
	"""
	Rscript db_build.R
	"""

	stub:
	"""
	mkdir -p orgfinder_db/fna
	touch orgfinder_db/fna/GCF_000000000.1.fna orgfinder_db/db.tsv
	"""
}


workflow ORGFINDER_DB_BUILD {
	main:
	def queries = channel.of(
		"taxon 'Pseudomonas aeruginosa'  --reference --assembly-level complete --include genome,seq-report",
		"taxon 'Acinetobacter baumannii' --reference --assembly-level complete --include genome,seq-report",
		"taxon 'Enterococcus'            --reference --assembly-level complete --include genome,seq-report",
		"taxon 'Staphylococcus'          --reference --assembly-level complete --include genome,seq-report",
		"taxon 'Streptococcus'           --reference --assembly-level complete --include genome,seq-report",
		"taxon 'Enterobacterales'        --reference --assembly-level complete --include genome,seq-report",
		"taxon 'Aeromonas'               --reference --assembly-level complete --include genome,seq-report",
		"taxon 'Myroides'                --reference --assembly-level complete --include genome,seq-report",
		"taxon 'Enterococcus faecalis'   --reference --include genome,seq-report",
		"taxon 'Citrobacter murliniae'   --reference --include genome,seq-report",
		'accession GCA_040096145.1', //E. intestinihominis	3133180	GCA_040096145.1	Genome of type strain CLA-AC-H004ᵀ
		'accession GCA_001875655.1', //E. hormaechei subsp. hormaechei	301105 / parent 158836	GCA_001875655.1	Genome of species type strain ATCC 49162ᵀ
		'accession GCA_001729745.1', //E. hormaechei subsp. hoffmannii	1812934	GCA_001729745.1	Complete genome of type strain DSM 14563ᵀ
		'accession GCF_048568405.1', //E. intestinihominis RefSeq representative	3133180	GCF_048568405.1	Our current reference
		'accession GCF_019048245.1'  //E. hormaechei RefSeq representative	158836	GCF_019048245.1	Our current reference
	)
	def datasets = NCBI_DATASET_DOWNLOAD_GENOME(queries.map({ q -> record(args: q) }))
	def ncbi_db = STEP2_MOLECULE_TYPES_TSV(STEP1_ACCESSIONS_TSV(datasets.collect()))

	def taxdump = UNTAR(channel.value(record(
		archive: file("https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"),
		dir_name: 'taxdump'
	)))

	def db_built = STEP3_DB_BUILD(ncbi_db, taxdump, file("${moduleDir}/assets/db_build.R"))

	emit:
	db: Value<OrgfinderDb> = db_built
}
