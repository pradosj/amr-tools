nextflow.enable.types = true

include { NANOPLOT           } from './modules/nanoplot'
include { FASTP              } from './modules/fastp'
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

// A file to put in the MultiQC folder under the given name
record NamedFile {
	name: String
	path: Path
}

// Gather files into a folder, each one under its given name
// (MultiQC derives sample names from the file names)
process MULTIQC_FOLDER {
    executor 'local'
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
