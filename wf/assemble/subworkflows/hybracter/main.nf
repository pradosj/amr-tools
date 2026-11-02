nextflow.enable.types = true

include { HYBRACTER as HYBRACTER_RUN } from './modules/hybracter'


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

// Assemble each sample with its long reads, plus its short reads when available
workflow HYBRACTER {
	take:
	short_reads: Channel<ShortReads>
	long_reads: Channel<LongReads>
	args: String?

	main:
	def asm_ch = HYBRACTER_RUN(
		short_reads
			.join(long_reads.map({ r -> record(meta: r.meta, long_reads: r.reads) }), by: 'meta', remainder: true)
			.map({ r -> r + record(args: args) })
	)

	emit:
	assemblies: Channel<Assembly> = asm_ch
}
