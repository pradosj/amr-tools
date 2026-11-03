nextflow.enable.types = true

// ---------------------------------------
//  Input reads and assemblies of the standalone pipelines, given as glob
//  patterns: the sample_id is taken from the file names
// ---------------------------------------

include { SAMTOOLS_FASTQ as ONT_BAM_TO_FASTQ } from './modules/samtools/fastq'

// ---------------------------------------

// Input reads given as glob patterns (params.reads)
record ReadsParams {
	ont: String?
	illumina: String?
}

// Input assemblies given as a glob pattern (params.assemblies)
record AssembliesParams {
	fasta: String?
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

// Assembly given as input: meta holds sample_id, and assembler_name 'none'
record Genome {
	meta: Map
	fasta: Path
	gff: Path?
	info: Path?
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

// ONT (one FASTQ, BAM or CRAM per sample, BAM/CRAM converted to FASTQ) and Illumina reads (single-end, or paired-end R1/R2)
workflow INPUT_READS {
	take:
	ont: String?
	illumina: String?

	main:
	def long_ch = channel.empty()
	if (ont) {
		def ont_ch = channel.fromPath(ont).map({ f ->
			record(meta: [sample_id: f.name.replaceAll(/\.(fastq\.gz|fq\.gz|bam|cram)$/, '')], reads: f)
		})
		def bam_re = '.*\\.(bam|cram)$'
		long_ch = ont_ch
			.filter({ r -> !(r.reads.name ==~ bam_re) })
			.mix(ONT_BAM_TO_FASTQ(
				ont_ch
					.filter({ r -> r.reads.name ==~ bam_re })
					.map({ r -> record(meta: r.meta, bam: r.reads, args: null) })
			))
	}
	def short_ch = channel.empty()
	if (illumina) {
		short_ch = channel.fromPath(illumina)
			.map({ f ->
				def id = (f.name ==~ ".*${illuminaRe()}") ? f.name.replaceAll(illuminaRe(), '') : f.name.replaceAll(/\.(fastq|fq)\.gz$/, '')
				tuple(id, f)
			})
			.groupBy()
			.map({ id, files -> shortReads(id, files) })
	}

	emit:
	short_reads: Channel<ShortReads> = short_ch
	long_reads: Channel<LongReads>   = long_ch
}

// Assemblies FASTA files
workflow INPUT_ASSEMBLIES {
	take:
	fasta: String?

	main:
	def asm_ch = channel.empty()
	if (fasta) {
		asm_ch = channel.fromPath(fasta).map({ f ->
			record(meta: [sample_id: f.name.replaceAll(/\.(fasta|fa|fna)$/, ''), assembler_name: 'none'], fasta: f)
		})
	}

	emit:
	assemblies: Channel<Genome> = asm_ch
}
