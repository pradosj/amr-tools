nextflow.enable.types = true

include { NANOPLOT                             } from './modules/nanoplot'
include { FASTP                                } from './modules/fastp'
include { MULTIQC                              } from './modules/multiqc'
include { MetaPath } from './subworkflows/input/main.nf'


// Illumina reads of a sample: read2 is null for single-end reads
record ShortReads {
	meta: Record?
	read1: Path
	read2: Path?
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


// QC of the reads of each sample, and a MultiQC report gathering all the samples
workflow READS_QC {
	take:
	short_reads: Channel<ShortReads>
	long_reads: Channel<MetaPath>     // one FASTQ per sample

	main:
	def nanoplot_ch = NANOPLOT(long_reads)
	def fastp_ch = FASTP(short_reads)

	// MultiQC
	def multiqc_files_ch = channel.empty()
		.mix(nanoplot_ch.map({ r -> record(name: "${r.meta.id}.nanostat", path: r.nanostat) })
		.mix(fastp_ch.json.map({ r -> record(name: "${r.meta.id}.fastp.json", path: r.path) }))
		.collect()
		.flatMap({ fs -> fs.isEmpty() ? [] : [fs.toList()] })  // no MultiQC when there are no reads
	def multiqc_ch = MULTIQC(
		MULTIQC_FOLDER(multiqc_files_ch)
		.map({d -> record(
		  meta: 'multiqc.html', 
		  db: d, 
		  config: file("${moduleDir}/assets/multiqc_config.yml")
		)})
	)

	emit:
	nanoplot: Channel<Record>    = nanoplot_ch
	fastp: Channel<Record>       = fastp_ch.outputs    // records with fields meta, html, json
	multiqc: Channel<Record>     = multiqc_ch    // records with fields meta, html
}
