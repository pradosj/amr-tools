nextflow.enable.types = true

include { ASSEMBLIES_QC; AssembliesQcParams } from './workflows/assemblies-qc/main.nf'
include { INPUT_READS; INPUT_ASSEMBLIES; ReadsParams; AssembliesParams    } from './subworkflows/input/main.nf'
include { DATABASES; DbParams; Db } from './subworkflows/db/main.nf'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'


// Input assemblies and reads given as glob patterns (the sample_id is taken from the file names), and ASSEMBLIES_QC options
params {
	db: DbParams = record(basedir: 'results/db', checkm2: 'checkm2', org_db: 'org_db', amrfinderplus: 'amrfinderplus')    // databases folder (built when missing), or tar archive of it (path or URL)
	assemblies: AssembliesParams = record(fasta: null)
	reads: ReadsParams = record(ont: null, illumina: null)
	assemblies_qc: AssembliesQcParams = record(checkm2: record(skip: false, args: ''))
}

workflow {
	main:
	// Validate parameters and print summary of supplied ones
	validateParameters()
	log.info(paramsSummaryLog(workflow))

	def assemblies = INPUT_ASSEMBLIES(params.assemblies?.fasta)
	def reads = INPUT_READS(params.reads?.ont, params.reads?.illumina)

	// Reads of the sample of each assembly
	def asm_keys = assemblies.map({ r -> record(meta: record(sample_id: r.meta.sample_id), asm_meta: r.meta) })
	def asm_short_reads = asm_keys.join(reads.short_reads, by: 'meta').map({ r -> record(meta: r.asm_meta, read1: r.read1, read2: r.read2) })
	def asm_long_reads = asm_keys.join(reads.long_reads, by: 'meta').map({ r -> record(meta: r.asm_meta, reads: r.reads) })

	def dbs = DATABASES(params.db, (params.assemblies_qc?.checkm2?.skip ? [] : [Db.CHECKM2]).toSet())
	def qc = ASSEMBLIES_QC(params.assemblies_qc, dbs.checkm2, assemblies, asm_short_reads, asm_long_reads)

	publish:
	db = dbs.built
	assemblies_checkm2 = qc.checkm2
	assemblies_multiqc = qc.multiqc
}

output {
	db: Channel<Record> {
		path { r -> r.dir >> "db/${r.name}" }    // databases built by this run: use --db=<output dir>/db in the next runs
	}
	assemblies_checkm2: Channel<Record> {
		path { r -> r.report >> "assemblies/${r.meta.sample_id}__${r.meta.assembler_name}.qc/checkm2_quality.tsv" }
	}
	assemblies_multiqc: Channel<Record> {
		path { r ->
			r.txt >> "reports/assemblies_multiqc.txt"
			r.xlsx >> "reports/assemblies_multiqc.xlsx"
		}
	}
}
