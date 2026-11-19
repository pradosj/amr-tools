nextflow.enable.types = true

include { GENOME_QC; GenomeQcParams } from './workflows/qc/main.nf'
include { INPUT_READS; INPUT_ASSEMBLIES; ReadsParams; AssembliesParams; MetaPath } from './subworkflows/input/main.nf'
include { DATABASES; DbParams; Db } from './subworkflows/db/main.nf'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'


// Input assemblies and reads given as glob patterns (the sample id is taken from the file names), and GENOME_QC options
params {
	db: DbParams = record(basedir: 'results/db', checkm2: 'checkm2', org_db: 'org_db', amrfinderplus: 'amrfinderplus')    // databases folder (built when missing), or tar archive of it (path or URL)
	assemblies: AssembliesParams = record(fasta: null)
	reads: ReadsParams = record(ont: null, illumina: null)
	genome_qc: GenomeQcParams = record(checkm2: record(skip: false, args: ''))
}

workflow {
	main:
	// Validate parameters and print summary of supplied ones
	validateParameters()
	log.info(paramsSummaryLog(workflow))

	def assemblies: Channel<MetaPath> = INPUT_ASSEMBLIES(params.assemblies?.fasta)
		.map({ r -> record(meta: record(id: r.meta.id, assembler_name: 'none'), path: r.path) })
	def reads = INPUT_READS(params.reads?.ont, params.reads?.illumina)

	// Each assembly with the reads of its sample
	def qc_input = assemblies
		.map({ r -> record(meta: record(id: r.meta.id), asm: r) })
		.join(reads.short_reads, by: 'meta', remainder: true)
		.join(reads.long_reads.map({ r -> record(meta: r.meta, reads: r.path) }), by: 'meta', remainder: true)
		.filter({ r -> r.asm != null })
		.map({ r -> record(
			meta: r.asm.meta,
			path: r.asm.path,
			short_reads: r.read1 ? record(meta: null, read1: r.read1, read2: r.read2) : null,
			long_reads: r.reads ? record(meta: null, reads: r.reads) : null
		) })

	def dbs = DATABASES(params.db, (params.genome_qc?.checkm2?.skip || params.genome_qc?.checkm2?.db ? [] : [Db.CHECKM2]).toSet())
	def qc = GENOME_QC(
		qc_input,
		dbs.checkm2.map({ db -> record(checkm2: record(
			skip: params.genome_qc?.checkm2?.skip,
			args: params.genome_qc?.checkm2?.args,
			db: params.genome_qc?.checkm2?.db ?: db    // the database given by --genome_qc.checkm2.db, else the one of the databases folder
		)) })
	)

	publish:
	db = dbs.built
	assemblies_checkm2 = qc.checkm2
	assemblies_qc_report = qc.report
	assemblies_multiqc = qc.multiqc
}

output {
	db: Channel<Record> {
		path { r -> r.dir >> "db/${r.name}" }    // databases built by this run: use --db=<output dir>/db in the next runs
	}
	assemblies_checkm2: Channel<Record> {
		path { r -> r.path >> "assemblies/${r.meta.id}__${r.meta.assembler_name}.qc/checkm2_quality.tsv" }
	}
	assemblies_qc_report: Channel<Record> {
		path { r -> r.html >> "assemblies/${r.meta.id}__${r.meta.assembler_name}.qc/qc_report.html" }
	}
	assemblies_multiqc: Channel<Record> {
		path { r ->
			r.txt >> "reports/assemblies_multiqc.txt"
			r.xlsx >> "reports/assemblies_multiqc.xlsx"
		}
	}
}
