nextflow.enable.types = true

include { BWA_INDEX   } from './modules/bwa/index'
include { BWA_MEM     } from './modules/bwa/mem'
include { PILON       } from './modules/pilon'


// Short reads sample: read2 is null for single-end reads
record ShortReads {
	meta: Map
	reads: Path
	read2: Path?
}

// Assembly of a sample: assembler output folder and assembly FASTA
record Assembly {
	meta: Map
	dir: Path
	fasta: Path
}

// Polish assemblies with short reads: map the short reads on the assembly, then run Pilon
workflow PILON_POLISH {
	take:
	assemblies: Channel<Assembly>
	short_reads: Channel<ShortReads>

	main:
	def fasta_ch = assemblies.map({ r -> record(meta: r.meta, fasta: r.fasta) })
	def index_ch = BWA_INDEX(fasta_ch)
	def bam_ch = BWA_MEM(index_ch.join(short_reads, by: 'meta'))
	def pilon_ch = PILON(fasta_ch.join(bam_ch, by: 'meta'))

	emit:
	polished: Channel<Assembly> = pilon_ch
}
