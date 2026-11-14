nextflow.enable.types = true

include { BWA_INDEX   } from './modules/bwa/index'
include { BWA_MEM     } from './modules/bwa/mem'
include { PILON       } from './modules/pilon'
include { MetaPath    } from './subworkflows/input/main.nf'


// Illumina reads of a sample: read2 is null for single-end reads
record ShortReads {
	meta: Record?
	read1: Path
	read2: Path?
}

// Assembly of a sample: assembler output folder and assembly FASTA
record Genome {
	meta: Record
	dir: Path
	fasta: Path
	info: Path?
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
	def pilon_ch = PILON(assemblies.map({ r -> record(meta: r.meta, fasta: r.path) }).join(bam_ch, by: 'meta'))

	emit:
	polished: Channel<Genome> = pilon_ch.outputs    // records with fields meta, dir, fasta
}
