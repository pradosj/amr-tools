nextflow.enable.types = true

include { READS_QC } from './workflows/reads-qc/main.nf'
include { READS_INPUT; ReadsParams                          } from './conf/inputs.nf'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'


// Input reads given as glob patterns, the sample_id is taken from the file names
params {
	reads: ReadsParams = record(ont: null, illumina: null)
}

workflow {
	main:
	// Validate parameters and print summary of supplied ones
	validateParameters()
	log.info(paramsSummaryLog(workflow))

	def reads = READS_INPUT(params.reads?.ont, params.reads?.illumina)
	def qc = READS_QC(reads.short_reads, reads.long_reads)

	publish:
	nanoplot = qc.nanoplot
	fastp = qc.fastp
	multiqc = qc.multiqc
}

output {
	nanoplot: Channel<Record> {
		path { r -> r.nanoplot >> "${r.meta.sample_id}/nanoplot" }
	}
	fastp: Channel<Record> {
		path { r ->
			r.json >> "${r.meta.sample_id}/fastp.json"
			r.html >> "${r.meta.sample_id}/fastp.html"
		}
	}
	multiqc: Channel<Record> {
		path { r -> r.html >> "multiqc.html" }
	}
}
