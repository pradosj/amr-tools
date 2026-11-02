nextflow.enable.types = true

include { NANOPLOT           } from './modules/nanoplot'
include { FASTP              } from './modules/fastp'
include { ORGANIZE_FILES     } from './modules/organize_files'
include { MULTIQC            } from './modules/multiqc'


workflow READS_QC {
	take:
	fqs_ch: Channel<Tuple>    // short reads: [meta, [R1, R2?]]
	fql_ch: Channel<Tuple>    // long reads:  [meta, fastq]

	main:
	def nanoplot_ch = NANOPLOT(fql_ch, null)
	def fastp_ch = FASTP(fqs_ch, null)

	// MultiQC
	def multiqc_files_ch = nanoplot_ch
		.map({ r -> tuple(r.nanostat, "${r.meta.sample_id}.nanostat") })
		.mix(fastp_ch.map({ r -> tuple(r.json, "${r.meta.sample_id}.fastp.json") }))
		.collect()
		.flatMap({ fs -> fs.isEmpty() ? [] : [tuple("multiqc.html", fs.toList())] })  // no MultiQC when there are no reads
	def multiqc_ch = MULTIQC(
		ORGANIZE_FILES(multiqc_files_ch),
		file("${moduleDir}/assets/multiqc_config.yml"),
		null
	)

	emit:
	nanoplot: Channel<Record>    = nanoplot_ch
	fastp: Channel<Record>       = fastp_ch
	multiqc_html: Channel<Path> = multiqc_ch.map({ r -> r.html })
}
