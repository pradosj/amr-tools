nextflow.enable.types = true

include { AMR_ANNOT; AmrParams } from './workflows/amr-annot/main.nf'
include { INPUT_READS; INPUT_ASSEMBLIES; ReadsParams; AssembliesParams } from './subworkflows/input/main.nf'
include { DATABASES; DbParams; Db } from './subworkflows/db/main.nf'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'


// Input assemblies and reads given as glob patterns (the sample_id is taken from the file names), and AMR_ANNOT options
params {
	db: DbParams = record(basedir: 'results/db', checkm2: 'checkm2', org_db: 'org_db', amrfinderplus: 'amrfinderplus')    // databases folder (built when missing), or tar archive of it (path or URL)
	assemblies: AssembliesParams = record(fasta: null)
	reads: ReadsParams = record(ont: null, illumina: null)
	amr: AmrParams = record(
		orgfinder:     record(skip: false, args: ''),
		speciator:     record(skip: true, args: ''),
		amrfinderplus: record(skip: false, args: ''),
		mobtyper:      record(skip: false, args: ''),
		plasmidfinder: record(skip: false, skip_longread: false, skip_shortread: false, args: ''),
		resfinder:     record(skip: false, skip_longread: false, skip_shortread: false, args: ''),
		cgemlst:       record(skip: false, args: ''),
		MLST:          record(skip: false, args: ''),    // args '' to autodetect the species
		prokka:        record(skip: true, args: '--kingdom Bacteria')
	)
}

workflow {
	main:
	// Validate parameters and print summary of supplied ones
	validateParameters()
	log.info(paramsSummaryLog(workflow))

	def assemblies = INPUT_ASSEMBLIES(params.assemblies?.fasta)
	def reads = INPUT_READS(params.reads?.ont, params.reads?.illumina)
	def dbs = DATABASES(params.db, ((params.amr?.orgfinder?.skip ? [] : [Db.ORG_DB]) + (params.amr?.amrfinderplus?.skip ? [] : [Db.AMRFINDERPLUS])).toSet())
	def amr_annot = AMR_ANNOT(params.amr, dbs.org_db, dbs.amrfinderplus, assemblies, reads.short_reads, reads.long_reads)

	publish:
	amr_multireport = amr_annot
	db = dbs.built
}

output {
	db: Channel<Record> {
		path { r -> r.dir >> "db/${r.name}" }    // databases built by this run: use --db=<output dir>/db in the next runs
	}
	amr_multireport: Channel<Record> {
		mode "copy"    // the folder holds symlinks to the task outputs
		path { r ->
			r.folder.resolve('amr') >> "amr"    // results of the AMR tools, as gathered for the multireport
			r.html >> "reports/report.html"
			r.xlsx >> "reports/multitable.xlsx"
			r.assemblies_tsv >> "reports/multitable_assemblies.tsv"
			r.contigs_tsv >> "reports/multitable_contigs.tsv"
			r.resistances_tsv >> "reports/multitable_resistances.tsv"
			r.plasmidfinder_hits_tsv >> "reports/multitable_plasmidfinder_hits.tsv"
		}
	}
}
