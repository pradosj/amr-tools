nextflow.enable.types = true

include { GENOME_QC; GenomeQcParams } from './workflows/qc/main.nf'
include { INPUT_READS; INPUT_GENOMES; ReadsParams; GenomesParams } from './subworkflows/input/main.nf'
include { DATABASES; DbParams; Db } from './subworkflows/db/main.nf'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'


// Input genomes and reads given as glob patterns (the sample id is taken from the file names), and GENOME_QC options
params {
	db: DbParams = record(basedir: 'results/db', checkm2: 'checkm2', org_db: 'org_db', amrfinderplus: 'amrfinderplus')    // databases folder (built when missing), or tar archive of it (path or URL)
	genomes: GenomesParams = record(fasta: null, json: null)
	reads: ReadsParams = record(ont: null, illumina: null)
	genome_qc: GenomeQcParams = record(checkm2: record(skip: false, args: ''))
}

workflow {
	main:
	// Validate parameters and print summary of supplied ones
	validateParameters()
	log.info(paramsSummaryLog(workflow))

	// the FASTA of each genome, with the asm.json of its assembler given with --genomes.json (null when not given)
	def genomes = INPUT_GENOMES(params.genomes?.fasta, params.genomes?.json)
		.map({ r -> record(meta: record(id: r.meta.id, assembler_name: 'none'), path: r.path, json: r.json) })
	def reads = INPUT_READS(params.reads?.ont, params.reads?.illumina)

	// Each genome with the reads of its sample
	def qc_input = genomes
		.map({ r -> record(meta: record(id: r.meta.id), asm: r) })
		.join(reads.short_reads, by: 'meta', remainder: true)
		.join(reads.long_reads.map({ r -> record(meta: r.meta, reads: r.path) }), by: 'meta', remainder: true)
		.filter({ r -> r.asm != null })
		.map({ r -> record(
			meta: r.asm.meta,
			path: r.asm.path,
			short_read1: r.read1,
			short_read2: r.read2,
			long_reads: r.reads,
			asm_info: r.asm.json
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
	genomes_checkm2 = qc.checkm2
	genomes_qc_report = qc.report
	genomes_multiqc = qc.multiqc
}

output {
	db: Channel<Record> {
		path { r -> r.dir >> "db/${r.name}" }    // databases built by this run: use --db=<output dir>/db in the next runs
	}
	genomes_checkm2: Channel<Record> {
		path { r -> r.path >> "genomes/${r.meta.id}__${r.meta.assembler_name}.qc/checkm2_quality.tsv" }
	}
	genomes_qc_report: Channel<Record> {
		path { r -> r.html >> "genomes/${r.meta.id}__${r.meta.assembler_name}.qc/qc_report.html" }
	}
	genomes_multiqc: Channel<Record> {
		path { r ->
			r.txt >> "reports/genomes_multiqc.txt"
			r.xlsx >> "reports/genomes_multiqc.xlsx"
		}
	}
}
