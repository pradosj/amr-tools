nextflow.enable.types = true

include { NANOPLOT                             } from './modules/nanoplot'
include { FASTP                                } from './modules/fastp'
include { MULTIQC                              } from './modules/multiqc'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

// Input reads given as glob patterns, the sample_id is taken from the file names
params {
	ont: String? = null
	illumina: String? = null
}

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

// A file to put in the MultiQC folder under the given name
record NamedFile {
	name: String
	path: Path
}

// Gather files into a folder, each one under its given name
// (MultiQC derives sample names from the file names)
process MULTIQC_FOLDER {
    label 'local'
    input:
    files: List<NamedFile>
    stage:
    stageAs(files.collect({ f -> f.path }), 'in/?/*')
    output:
    file('multiqc_folder')
    script:
    def moves = files.collect({ f -> "mv '${f.path}' 'multiqc_folder/${f.name}'" }).join('\n')
    """
    mkdir multiqc_folder
    ${moves}
    """
}


// Paired-end Illumina file names end with _1, _2, _R1 or _R2 (optionally followed by _001)
def illuminaRe() -> String {
	return '_(R?)([12])(_[0-9][0-9][0-9])?\\.(fq|fastq)\\.gz$'
}

// Short reads of a sample from its files: R1 and R2 when paired-end, a single file otherwise
def shortReads(id: String, files: Bag<Path>) -> ShortReads {
	def fs = files.toList().toSorted({ f -> f.name })
	if (fs.size() > 2) {
		def msg = "Sample ${id}: more than 2 short reads files: ${fs}"
		log.error(msg)    // errors raised in operators are only reported in .nextflow.log
		error(msg)
	}
	return record(meta: [sample_id: id], reads: fs[0], read2: fs.size() == 2 ? fs[1] : null)
}


// QC of the reads of each sample, and a MultiQC report gathering all the samples
workflow READS_QC {
	take:
	short_reads: Channel<ShortReads>
	long_reads: Channel<LongReads>

	main:
	def nanoplot_ch = NANOPLOT(long_reads)
	def fastp_ch = FASTP(short_reads)

	// MultiQC
	def multiqc_files_ch = nanoplot_ch
		.map({ r -> record(name: "${r.meta.sample_id}.nanostat", path: r.nanostat) })
		.mix(fastp_ch.map({ r -> record(name: "${r.meta.sample_id}.fastp.json", path: r.json) }))
		.collect()
		.flatMap({ fs -> fs.isEmpty() ? [] : [fs.toList()] })  // no MultiQC when there are no reads
	def multiqc_config = file("${moduleDir}/assets/multiqc_config.yml")
	def multiqc_ch = MULTIQC(
		MULTIQC_FOLDER(multiqc_files_ch).map({ d -> record(meta: 'multiqc.html', db: d, config: multiqc_config) })
	)

	emit:
	nanoplot: Channel<Record>    = nanoplot_ch
	fastp: Channel<Record>       = fastp_ch
	multiqc: Channel<Record>     = multiqc_ch    // records with fields meta, html
}


workflow {
	main:
	// Validate parameters and print summary of supplied ones
	validateParameters()
	log.info(paramsSummaryLog(workflow))

	def long_ch = channel.empty()
	if (params.ont) {
		long_ch = channel.fromPath(params.ont).map({ f ->
			record(meta: [sample_id: f.name.replaceAll(/\.(fastq|fq)\.gz$/, '')], reads: f)
		})
	}
	def short_ch = channel.empty()
	if (params.illumina) {
		short_ch = channel.fromPath(params.illumina)
			.map({ f ->
				def id = (f.name ==~ ".*${illuminaRe()}") ? f.name.replaceAll(illuminaRe(), '') : f.name.replaceAll(/\.(fastq|fq)\.gz$/, '')
				tuple(id, f)
			})
			.groupBy()
			.map({ id, files -> shortReads(id, files) })
	}

	def qc = READS_QC(short_ch, long_ch)

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
