nextflow.enable.types = true

// ---------------------------------------
//  Input reads and assemblies of the standalone pipelines, given as comma-separated lists of glob
//  patterns and URLs: the sample_id is taken from the file names
// ---------------------------------------

include { SAMTOOLS_FASTQ as ONT_BAM_TO_FASTQ } from './modules/samtools/fastq'
include { GZIP_DECOMPRESS                    } from './modules/gzip'

// ---------------------------------------

// Input reads given as comma-separated lists of glob patterns and URLs (params.reads)
record ReadsParams {
	ont: String?
	illumina: String?
}

// Input assemblies given as a comma-separated list of glob patterns and URLs (params.assemblies)
record AssembliesParams {
	fasta: String?
}

// Illumina reads of a sample: read2 is null for single-end reads
record ShortReads {
	meta: Record?
	read1: Path
	read2: Path?
}

// ONT reads of a sample: one FASTQ
record LongReads {
	meta: Record?
	reads: Path
}

// Sequencing reads of a sample: high quality ONT reads and/or Illumina reads
record SeqReads {
	meta: Record?
	ont_hq: LongReads?
	illumina: ShortReads?
}

// Assembly given as input: meta holds sample_id, and assembler_name 'none'
record Genome {
	meta: Record
	fasta: Path
	gff: Path?
	info: Path?
}

// A file and its meta (generic: e.g. the output of a module)
record MetaPath {
	meta: Record?
	path: Path
}

// A meta holding the (non null) id of its file, other fields allowed
record IdMeta {
	id: String
}

// A file and its meta holding its id
record IdPath {
	meta: IdMeta
	path: Path
}

// The record r ({ meta, path }, e.g. a module output, whose meta is typed Record) as an IdPath: its meta reduced to its
// id (error when the meta or its id is null)
def asIdPath(r: Record) -> IdPath {
	if (r.meta?.id == null)
		error("asIdPath: the meta of ${r.path} has no id")
	return record(meta: record(id: r.meta.id as String), path: r.path as Path)
}

// Files of a comma-separated list of glob patterns and URLs, each pattern sorted (commas inside the braces of a glob do not separate).
// Remote files are not listed: a URL must name a single file, without glob characters.
def pathList(patterns: String) -> List<Path> {
	def ps = patternList(patterns)
	def bad = ps.findAll({ p -> p ==~ /^[a-z]+:\/\/.*[*{\[].*/ })    // ? is the start of a query string
	if (bad) {
		error("Glob patterns are not supported in URLs, list each file: ${bad.join(', ')}")
	}
	return ps.collectMany({ p -> files(p).toSorted() }).toList()
}

// The glob patterns and URLs of a comma-separated list (commas inside the braces of a glob do not separate)
def patternList(patterns: String) -> List<String> {
	return patterns.replaceAll(/,(?![^{}]*\})/, "\n").tokenize("\n").collect({ p -> p.trim() }).toList()
}

// A sequence of a FASTA file to split (its FASTA entry) written to its own FASTA file
process SPLIT_MULTIFASTA {
	input:
	record(meta: IdMeta, text: String)
	output:
	record(meta: meta, path: file('sequence.fasta'))
	exec:    // native task: the file is written by Nextflow, without a shell script
	task.workDir.resolve('sequence.fasta').text = text
}

// Sample id of a FASTA file: its name without the .fasta, .fa or .fna extension, optionally followed by .gz
def fastaSampleId(name: String) -> String {
	return name.replaceAll(/\.(fasta|fa|fna)(\.gz)?$/, '')
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
	return record(meta: record(sample_id: id), read1: fs[0], read2: fs.size() == 2 ? fs[1] : null)
}

// ONT (one FASTQ, BAM or CRAM per sample, BAM/CRAM converted to FASTQ) and Illumina reads (single-end, or paired-end R1/R2):
// by type of reads (short_reads, long_reads), and both types together by sample (seq_reads)
workflow INPUT_READS {
	take:
	ont: String?
	illumina: String?

	main:
	def long_ch = channel.empty()
	if (ont) {
		def ont_ch: Channel<LongReads> = channel.fromList(pathList(ont)).map({ f ->
			record(meta: record(sample_id: f.name.replaceAll(/\.(fastq\.gz|fq\.gz|bam|cram)$/, '')), reads: f)
		})
		def bam_re = '.*\\.(bam|cram)$'
		long_ch = ont_ch
			.filter({ r -> !(r.reads.name ==~ bam_re) })
			.mix(ONT_BAM_TO_FASTQ(
				ont_ch
					.filter({ r -> r.reads.name ==~ bam_re })
					.map({ r -> record(meta: r.meta, bam: r.reads) })
			))
	}
	def short_ch = channel.empty()
	if (illumina) {
		short_ch = channel.fromList(pathList(illumina))
			.map({ f ->
				def id = (f.name ==~ ".*${illuminaRe()}") ? f.name.replaceAll(illuminaRe(), '') : f.name.replaceAll(/\.(fastq|fq)\.gz$/, '')
				tuple(id, f)
			})
			.groupBy()
			.map({ id, files -> shortReads(id, files) })
	}

	def seq_ch = short_ch
		.map({ r -> record(meta: r.meta, illumina: record(read1: r.read1, read2: r.read2)) })
		.join(long_ch.map({ r -> record(meta: r.meta, ont_hq: record(reads: r.reads)) }), by: 'meta', remainder: true)
		.map({ r -> record(meta: r.meta, ont_hq: r.ont_hq, illumina: r.illumina) })

	emit:
	seq_reads: Channel<SeqReads>     = seq_ch
	short_reads: Channel<ShortReads> = short_ch
	long_reads: Channel<LongReads>   = long_ch
}

// FASTA files of a comma-separated list of paths, glob patterns and URLs: gzipped files are decompressed.
// IdPath records: meta is record(id), the id being the sample id of the file (see fastaSampleId()).
// A pattern suffixed by `!` (e.g. `plasmids.fasta!`) gives files to split: each of their sequences is a genome, written
// to its own FASTA file, with id `<sample id of the file>:<sequence id>`.
workflow INPUT_FASTA {
	take:
	fasta: String?

	main:
	def fasta_ch = channel.empty()
	if (fasta) {
		def ps = patternList(fasta)
		def split_ch = SPLIT_MULTIFASTA(
			channel.fromList(pathList(ps.findAll({ p -> p.endsWith('!') }).collect({ p -> p.substring(0, p.length() - 1) }).join(',')))
				.flatMap({ f -> f.splitFasta(record: [id: true, text: true]).collect({ s ->     // .gz files are read decompressed
					record(meta: record(id: "${fastaSampleId(f.name)}:${s.id}"), text: s.text)
				}) })
		)
		def in_ch: Channel<IdPath> = channel.fromList(pathList(ps.findAll({ p -> !p.endsWith('!') }).join(',')))
			.map({ f -> record(meta: record(id: fastaSampleId(f.name)), path: f) })
		fasta_ch = in_ch
			.filter({ r -> !r.path.name.endsWith('.gz') })
			.mix(GZIP_DECOMPRESS(
				in_ch
					.filter({ r -> r.path.name.endsWith('.gz') })
			).map({ r -> asIdPath(record(meta: r.meta, path: r.file)) }))
			.mix(split_ch)
	}

	emit:
	fasta_files: Channel<IdPath> = fasta_ch
}

// Assemblies FASTA files
workflow INPUT_ASSEMBLIES {
	take:
	fasta: String?

	main:
	def asm_ch = INPUT_FASTA(fasta).map({ r -> record(meta: record(sample_id: r.meta.id, assembler_name: 'none'), fasta: r.path) })

	emit:
	assemblies: Channel<Genome> = asm_ch
}
