nextflow.enable.types = true

include { AMR_ANNOT; AmrParams } from './workflows/amr-annot/main.nf'
include { READS_INPUT; ASSEMBLIES_INPUT; ReadsParams; AssembliesParams } from './conf/inputs.nf'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'


// Input assemblies and reads given as glob patterns (the sample_id is taken from the file names), and AMR_ANNOT options
params {
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

	def assemblies = ASSEMBLIES_INPUT(params.assemblies?.fasta)
	def reads = READS_INPUT(params.reads?.ont, params.reads?.illumina)
	def report = AMR_ANNOT(params.amr, assemblies, reads.short_reads, reads.long_reads)

	publish:
	amr_multireport = report
}

output {
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
