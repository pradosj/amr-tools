
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
		def nanoplot_ch = NANOPLOT(fql_ch, null)

		// Run fastp once on each FASTQ-pair, and then expand to inputs sharing the same FASTQ
		def fastp_ch = FASTP(fqs_ch, null)

		// MultiQC
		ORGANIZE_FILES(
			channel.empty().mix(
				nanoplot_ch.map({r -> [r.nanostat,"${r.meta.sample_id}.nanostat"]}),
				fastp_ch.map({r -> [r.json,"${r.meta.sample_id}.fastp.json"]})
			)
			.collect({f -> [f]})
			.map({fs -> ["multiqc.html",fs]})
		)
		def multiqc_ch = MULTIQC(ORGANIZE_FILES.out, file("${moduleDir}/assets/multiqc_config.yml"), null)
	emit:
		long_nanoplot     = nanoplot_ch.map({r -> [r.meta, r.nanoplot]})
		long_nanostat     = nanoplot_ch.map({r -> [r.meta, r.nanostat]})
		short_fastp_json  = fastp_ch.map({r -> [r.meta, r.json]})
		short_fastp_html  = fastp_ch.map({r -> [r.meta, r.html]})
		multiqc_html      = multiqc_ch.map({r -> r.html})
}
