nextflow.enable.types = true

include { FLYE                          } from './modules/flye'
include { LONG_FLYE_MEDAKA              } from './subworkflows/flye_medaka_pilon'
include { HYBRID_FLYE_MEDAKA_PILON      } from './subworkflows/flye_medaka_pilon'
include { UNICYCLER                     } from './modules/unicycler'
include { SPADES                        } from './modules/spades'
include { HYBRACTER                     } from './modules/hybracter'
include { MetaPath                      } from './subworkflows/input/main.nf'


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

// Outputs of an assembler: the assembly FASTA and the folder of the assembler outputs
record AssemblyOutputs {
	meta: Record?
	fasta: Path
	dir: Path
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

	def fasta_ch: Channel<MetaPath> = channel.empty()
	def outputs_ch: Channel<AssemblyOutputs> = channel.empty()
	if (name == 'long_flye') {
		def asm = FLYE(long_input.map({ r -> record(meta: r.meta, path: r.long_reads, args: opts.args?.flye) }))
		fasta_ch = asm.fasta
		outputs_ch = asm.outputs
	} else if (name == 'long_flye_medaka') {
		def asm = LONG_FLYE_MEDAKA(long_input.map({ r -> r + record(args: opts.args) }))
		fasta_ch = asm.fasta
		outputs_ch = asm.outputs
	} else if (name == 'long_unicycler') {
		def asm = UNICYCLER(long_input.map({ r -> r + record(args: opts.args?.unicycler) }))
		fasta_ch = asm.fasta
		outputs_ch = asm.outputs
	} else if (name == 'long_hybracter') {
		def asm = HYBRACTER(long_input.map({ r -> r + record(args: opts.args?.hybracter) }))
		fasta_ch = asm.fasta
		outputs_ch = asm.outputs
	} else if (name == 'short_unicycler') {
		def asm = UNICYCLER(short_input.map({ r -> r + record(args: opts.args?.unicycler) }))
		fasta_ch = asm.fasta
		outputs_ch = asm.outputs
	} else if (name == 'short_spades') {
		def asm = SPADES(short_input.map({ r -> r + record(args: opts.args?.spades) }))
		fasta_ch = asm.fasta
		outputs_ch = asm.outputs
	} else if (name == 'hybrid_unicycler') {
		// samples with a single type of reads are assembled with it alone
		def asm = UNICYCLER(hybrid_input.map({ r -> r + record(args: opts.args?.unicycler) }))
		fasta_ch = asm.fasta
		outputs_ch = asm.outputs
	} else if (name == 'hybrid_hybracter') {
		// samples with long reads only are assembled with `hybracter long-single`,
		// samples with short reads only are dropped (hybracter needs long reads)
		def asm = HYBRACTER(hybrid_input.filter({ r -> r.long_reads != null }).map({ r -> r + record(args: opts.args?.hybracter) }))
		fasta_ch = asm.fasta
		outputs_ch = asm.outputs
	} else if (name == 'hybrid_flye_medaka_pilon') {
		// samples with short reads only are dropped (flye needs long reads)
		def asm = HYBRID_FLYE_MEDAKA_PILON(hybrid_input.filter({ r -> r.long_reads != null }).map({ r -> r + record(args: opts.args) }))
		fasta_ch = asm.fasta
		outputs_ch = asm.outputs
	} else {
		error("Unknown assembler name: ${name}")
	}

	emit:
	// meta holds id (sample id) and assembler_name
	fasta: Channel<MetaPath> = fasta_ch.map({ r -> record(meta: r.meta + record(assembler_name: name), path: r.path) })    // the assembly FASTA
	outputs: Channel<AssemblyOutputs> = outputs_ch.map({ r -> r + record(meta: r.meta + record(assembler_name: name)) })    // records with fields meta, fasta, dir (assembler outputs)
}
