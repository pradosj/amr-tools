nextflow.enable.types = true

include { NANOPLOT           } from './modules/nanoplot'
include { FASTP              } from './modules/fastp'
include { ORGANIZE_FILES     } from './modules/organize_files'
include { MULTIQC            } from './modules/multiqc'


// Short reads sample: read2 is null for single-end reads
record ShortReads {
	meta: Map
	reads: Path
	read2: Path?
}

// Long reads sample: one FASTQ per sample
record LongReads {
	meta: Map
	reads: Path
}

workflow READS_QC {
	take:
	short_reads: Channel<ShortReads>
	long_reads: Channel<LongReads>

	main:
	def nanoplot_ch = NANOPLOT(long_reads)
	def fastp_ch = FASTP(short_reads)

	// MultiQC
	def multiqc_files_ch = nanoplot_ch
		.map({ r -> tuple(r.nanostat, "${r.meta.sample_id}.nanostat") })
		.mix(fastp_ch.map({ r -> tuple(r.json, "${r.meta.sample_id}.fastp.json") }))
		.collect()
		.flatMap({ fs -> fs.isEmpty() ? [] : [tuple("multiqc.html", fs.toList())] })  // no MultiQC when there are no reads
	def multiqc_config = file("${moduleDir}/assets/multiqc_config.yml")
	def multiqc_ch = MULTIQC(
		ORGANIZE_FILES(multiqc_files_ch).map({ t -> record(meta: t[0], db: t[1], config: multiqc_config) })
	)

	emit:
	nanoplot: Channel<Record>    = nanoplot_ch
	fastp: Channel<Record>       = fastp_ch
	multiqc_html: Channel<Path> = multiqc_ch.map({ r -> r.html })
}
