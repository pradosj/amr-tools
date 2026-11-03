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


workflow ASSEMBLE {
	take:
	opts: AssemblerParams
	short_reads: Channel<ShortReads>
	long_reads: Channel<LongReads>

	main:
	def name = opts.name ?: 'long_flye_medaka'
	def args = opts.args
	def asm_ch = channel.empty()
	if (name == 'long_flye') {
		asm_ch = LONG_FLYE(long_reads, args?.flye)
	} else if (name == 'long_flye_medaka') {
		asm_ch = LONG_FLYE_MEDAKA(long_reads, args?.flye, args?.medaka)
	} else if (name == 'long_unicycler') {
		asm_ch = UNICYCLER(long_reads.map({ r -> record(meta: r.meta, reads: null, read2: null, long_reads: r.reads, args: args?.unicycler) }))
	} else if (name == 'long_hybracter') {
		asm_ch = HYBRACTER(long_reads.map({ r -> record(meta: r.meta, reads: null, read2: null, long_reads: r.reads, args: args?.hybracter) }))
	} else if (name == 'short_unicycler') {
		asm_ch = UNICYCLER(short_reads.map({ r -> r + record(long_reads: null, args: args?.unicycler) }))
	} else if (name == 'short_spades') {
		asm_ch = SPADES(short_reads.map({ r -> r + record(long_reads: null, args: args?.spades) }))
	} else if (name == 'hybrid_unicycler') {
		// samples with a single type of reads are assembled with it alone
		asm_ch = UNICYCLER(
			short_reads
				.join(long_reads.map({ r -> record(meta: r.meta, long_reads: r.reads) }), by: 'meta', remainder: true)
				.map({ r -> r + record(args: args?.unicycler) })
		)
	} else if (name == 'hybrid_hybracter') {
		// samples with long reads only are assembled with `hybracter long-single`,
		// samples with short reads only are dropped (hybracter needs long reads)
		asm_ch = HYBRACTER(
			short_reads
				.join(long_reads.map({ r -> record(meta: r.meta, long_reads: r.reads) }), by: 'meta', remainder: true)
				.filter({ r -> r.long_reads != null })
				.map({ r -> r + record(args: args?.hybracter) })
		)
	} else if (name == 'hybrid_flye_medaka_pilon') {
		asm_ch = HYBRID_FLYE_MEDAKA_PILON(short_reads, long_reads, args?.flye, args?.medaka)
	} else {
		error("Unknown assembler name: ${name}")
	}

	emit:
	assemblies: Channel<Assembly> = asm_ch    // meta holds sample_id and assembler_name
		.map({ r -> r + record(meta: r.meta + [assembler_name: name]) })
}
