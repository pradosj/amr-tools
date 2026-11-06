nextflow.enable.types = true

include { LONG_FLYE                     } from './subworkflows/flye_medaka_pilon'
include { LONG_FLYE_MEDAKA              } from './subworkflows/flye_medaka_pilon'
include { HYBRID_FLYE_MEDAKA_PILON      } from './subworkflows/flye_medaka_pilon'
include { UNICYCLER                     } from './modules/unicycler'
include { SPADES                        } from './modules/spades'
include { HYBRACTER                     } from './modules/hybracter'


// Arguments given to each assembly tool. Fields are nullable: modules fall back to their own defaults.
record AssemblerArgs {
	flye: String?
	medaka: String?
	spades: String?
	unicycler: String?
	hybracter: String?
}

// Options of ASSEMBLE (params.assembler)
record AssemblerParams {
	name: String?              // assembly strategy (default: long_flye_medaka)
	args: AssemblerArgs?
}

// Input of an assembler: the reads of a sample (short_read2 is null for single-end reads) and the tools arguments
record AssemblyInput {
	meta: Record?
	short_read1: Path?
	short_read2: Path?
	long_reads: Path?
	args: AssemblerArgs?
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

// Assembly of a sample: assembler output folder and assembly FASTA
record Genome {
	meta: Record
	dir: Path
	fasta: Path
	info: Path?
}


workflow ASSEMBLE {
	take:
	opts: AssemblerParams
	short_reads: Channel<ShortReads>
	long_reads: Channel<LongReads>

	main:
	def name = opts.name ?: 'long_flye_medaka'

	// Inputs of the assemblers: long reads only, short reads only, or both (hybrid), the arguments are added for each assembler
	def long_input = long_reads.map({ r -> record(meta: r.meta, long_reads: r.reads) })
	def short_input = short_reads.map({ r -> record(meta: r.meta, short_read1: r.read1, short_read2: r.read2) })
	def hybrid_input = short_reads
		.join(long_reads.map({ r -> record(meta: r.meta, long_reads: r.reads) }), by: 'meta', remainder: true)
		.map({ r -> record(meta: r.meta, short_read1: r.read1, short_read2: r.read2, long_reads: r.long_reads) })

	def asm_ch = channel.empty()
	if (name == 'long_flye') {
		asm_ch = LONG_FLYE(long_input.map({ r -> r + record(args: opts.args) }))
	} else if (name == 'long_flye_medaka') {
		asm_ch = LONG_FLYE_MEDAKA(long_input.map({ r -> r + record(args: opts.args) }))
	} else if (name == 'long_unicycler') {
		asm_ch = UNICYCLER(long_input.map({ r -> r + record(args: opts.args?.unicycler) }))
	} else if (name == 'long_hybracter') {
		asm_ch = HYBRACTER(long_input.map({ r -> r + record(args: opts.args?.hybracter) }))
	} else if (name == 'short_unicycler') {
		asm_ch = UNICYCLER(short_input.map({ r -> r + record(args: opts.args?.unicycler) }))
	} else if (name == 'short_spades') {
		asm_ch = SPADES(short_input.map({ r -> r + record(args: opts.args?.spades) }))
	} else if (name == 'hybrid_unicycler') {
		// samples with a single type of reads are assembled with it alone
		asm_ch = UNICYCLER(hybrid_input.map({ r -> r + record(args: opts.args?.unicycler) }))
	} else if (name == 'hybrid_hybracter') {
		// samples with long reads only are assembled with `hybracter long-single`,
		// samples with short reads only are dropped (hybracter needs long reads)
		asm_ch = HYBRACTER(hybrid_input.filter({ r -> r.long_reads != null }).map({ r -> r + record(args: opts.args?.hybracter) }))
	} else if (name == 'hybrid_flye_medaka_pilon') {
		// samples with short reads only are dropped (flye needs long reads)
		asm_ch = HYBRID_FLYE_MEDAKA_PILON(hybrid_input.filter({ r -> r.long_reads != null }).map({ r -> r + record(args: opts.args) }))
	} else {
		error("Unknown assembler name: ${name}")
	}

	emit:
	assemblies: Channel<Genome> = asm_ch    // meta holds sample_id and assembler_name
		.map({ r -> r + record(meta: r.meta + record(assembler_name: name)) })
}
