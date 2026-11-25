nextflow.enable.types = true

include { BWA_INDEX   } from './modules/aln/bwa/index'
include { BWA_MEM     } from './modules/aln/bwa/mem'
include { PILON       } from './modules/asm/pilon'
include { MetaPath    } from './subworkflows/input/main.nf'


// Illumina reads of a sample: read2 is null for single-end reads
record ShortReads {
	meta: Record?
	read1: Path
	read2: Path?
}

// An assembly: its FASTA and the folder of the tool outputs
record Assembly {
	meta: Record?
	path: Path
	dir: Path
}

// Polish assemblies with short reads: map the short reads on the assembly, then run Pilon
workflow PILON_POLISH {
	take:
	assemblies: Channel<MetaPath>    // the FASTA of each assembly
	short_reads: Channel<ShortReads>

	main:
	def index_ch = BWA_INDEX(assemblies)
		.map({ i -> record(meta: i.meta, ref: i) })
	def bam_ch = BWA_MEM(short_reads.join(index_ch, by: 'meta')).outputs
	def pilon_ch = PILON(assemblies.join(bam_ch, by: 'meta'))

	emit:
	polished: Channel<Assembly> = pilon_ch    // path: pilon.fasta, dir: the Pilon folder; a single emit: the workflow returns this channel
}
