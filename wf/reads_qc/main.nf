
include { NANOPLOT           } from './modules/nanoplot'
include { FASTP              } from './modules/fastp'
include { ORGANIZE_FILES     } from './modules/organize_files'
include { MULTIQC            } from './modules/multiqc'


workflow READS_QC {
	take:
		fqs_ch    // channel: [ val(meta), path(short_reads) ]	
		fql_ch    // channel: [ val(meta), path(long_reads) ]
	main:
		// Run nanoplot once on each FASTQ, and then expand to inputs sharing the same FASTQ
		NANOPLOT(fql_ch)

		// Run fastp once on each FASTQ-pair, and then expand to inputs sharing the same FASTQ
		def fastp_ch = FASTP(fqs_ch, null)

		// MultiQC
		ORGANIZE_FILES(
			channel.empty().mix(
				NANOPLOT.out.nanostat.map({meta,file -> [file,"${meta.sample_id}.nanostat"]}),
				fastp_ch.map({r -> [r.json,"${r.meta.sample_id}.fastp.json"]})
			)
			.collect({f -> [f]})
			.map({fs -> ["multiqc.html",fs]})
		)
		MULTIQC(ORGANIZE_FILES.out,channel.value(file("${moduleDir}/assets/multiqc_config.yml")))
	emit:
		long_nanoplot     = NANOPLOT.out.nanoplot
		long_nanostat     = NANOPLOT.out.nanostat
		short_fastp_json  = fastp_ch.map({r -> [r.meta, r.json]})
		short_fastp_html  = fastp_ch.map({r -> [r.meta, r.html]})
		multiqc_html      = MULTIQC.out.html.map({_m,x -> x})
}
