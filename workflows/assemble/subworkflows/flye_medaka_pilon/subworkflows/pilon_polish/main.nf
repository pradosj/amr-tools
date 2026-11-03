nextflow.enable.types = true

include { BWA_INDEX   } from './modules/bwa/index'
include { BWA_MEM     } from './modules/bwa/mem'
include { PILON       } from './modules/pilon'


// Illumina reads of a sample: read2 is null for single-end reads
record ShortReads {
	meta: Map?
	read1: Path
	read2: Path?
}

// Assembly of a sample: assembler output folder and assembly FASTA
record Genome {
	meta: Map
	dir: Path
	fasta: Path
	info: Path?
}

// Polish assemblies with short reads: map the short reads on the assembly, then run Pilon
workflow PILON_POLISH {
	take:
	assemblies: Channel<Genome>
	short_reads: Channel<ShortReads>

	main:
	def fasta_ch = assemblies.map({ r -> record(meta: r.meta, fasta: r.fasta) })
	def index_ch = BWA_INDEX(fasta_ch)
	def bam_ch = BWA_MEM(index_ch.join(short_reads, by: 'meta'))
	def pilon_ch = PILON(fasta_ch.join(bam_ch, by: 'meta'))

	emit:
	polished: Channel<Genome> = pilon_ch
}
