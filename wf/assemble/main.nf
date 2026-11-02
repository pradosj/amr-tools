nextflow.enable.types = true

include { LONG_FLYE                     } from './subworkflows/flye_medaka_pilon'
include { LONG_FLYE_MEDAKA              } from './subworkflows/flye_medaka_pilon'
include { UNICYCLER as LONG_UNICYCLER   } from './subworkflows/unicycler'
include { HYBRACTER as LONG_HYBRACTER   } from './subworkflows/hybracter'
include { UNICYCLER as SHORT_UNICYCLER  } from './subworkflows/unicycler'
include { SPADES    as SHORT_SPADES     } from './subworkflows/spades'
include { UNICYCLER as HYBRID_UNICYCLER } from './subworkflows/unicycler'
include { HYBRACTER as HYBRID_HYBRACTER } from './subworkflows/hybracter'
include { HYBRID_FLYE_MEDAKA_PILON      } from './subworkflows/flye_medaka_pilon'


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
	name: String               // assembly strategy, e.g. long_flye_medaka
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
	def args = opts.args
	def asm_ch = channel.empty()
	if (opts.name == 'long_flye') {
		asm_ch = LONG_FLYE(long_reads, args?.flye)
	} else if (opts.name == 'long_flye_medaka') {
		asm_ch = LONG_FLYE_MEDAKA(long_reads, args?.flye, args?.medaka)
	} else if (opts.name == 'long_unicycler') {
		asm_ch = LONG_UNICYCLER(channel.empty(), long_reads, args?.unicycler)
	} else if (opts.name == 'long_hybracter') {
		asm_ch = LONG_HYBRACTER(channel.empty(), long_reads, args?.hybracter)
	} else if (opts.name == 'short_unicycler') {
		asm_ch = SHORT_UNICYCLER(short_reads, channel.empty(), args?.unicycler)
	} else if (opts.name == 'short_spades') {
		asm_ch = SHORT_SPADES(short_reads, channel.empty(), args?.spades)
	} else if (opts.name == 'hybrid_unicycler') {
		asm_ch = HYBRID_UNICYCLER(short_reads, long_reads, args?.unicycler)
	} else if (opts.name == 'hybrid_hybracter') {
		asm_ch = HYBRID_HYBRACTER(short_reads, long_reads, args?.hybracter)
	} else if (opts.name == 'hybrid_flye_medaka_pilon') {
		asm_ch = HYBRID_FLYE_MEDAKA_PILON(short_reads, long_reads, args?.flye, args?.medaka)
	} else {
		error("Unknown assembler name: ${opts.name}")
	}

	emit:
	assemblies: Channel<Assembly> = asm_ch
}
