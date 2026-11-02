nextflow.enable.types = true

include { SPADES as SPADES_RUN } from './modules/spades'


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

// Assembly of a sample: assembler output folder and assembly FASTA
record Assembly {
	meta: Map
	dir: Path
	fasta: Path
}

// Assemble each sample with its short reads and/or long reads
workflow SPADES {
	take:
	short_reads: Channel<ShortReads>
	long_reads: Channel<LongReads>
	args: String?

	main:
	def asm_ch = SPADES_RUN(
		short_reads
			.join(long_reads.map({ r -> record(meta: r.meta, long_reads: r.reads) }), by: 'meta', remainder: true)
			.map({ r -> r + record(args: args) })
	)

	emit:
	assemblies: Channel<Assembly> = asm_ch
}
