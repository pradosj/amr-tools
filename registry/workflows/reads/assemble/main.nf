nextflow.enable.types = true

include { FLYE                          } from './modules/asm/flye'
include { LONG_FLYE_MEDAKA              } from './subworkflows/flye_medaka_pilon'
include { HYBRID_FLYE_MEDAKA_PILON      } from './subworkflows/flye_medaka_pilon'
include { UNICYCLER                     } from './modules/asm/unicycler'
include { SPADES                        } from './modules/asm/spades'
include { HYBRACTER                     } from './modules/asm/hybracter'


// Arguments given to each assembly tool. Fields are nullable: modules fall back to their own defaults.
record AssemblerArgs {
	flye: String?
	medaka: String?
	spades: String?
	unicycler: String?
	hybracter: String?
}

// Options of the assemblies (params.assembler of the pipelines, turned into AssemblyJob records by them)
record AssemblerParams {
	name: String?              // assembly strategies, comma-separated (default: long_flye_medaka): each sample is assembled by each of them
	args: AssemblerArgs?
}


// An assembly: its FASTA and the folder of the assembler outputs (the output of each assembler)
record Assembly {
	meta: Record?
	path: Path
	dir: Path
}

// An assembly to do: the reads of a sample (null when missing), the strategy and the tools arguments
record AssemblyJob {
	meta: Record           // id (sample id) and assembler_name (names the assembly in the outputs, unique to each job)
	assembler: String      // assembly strategy: long_flye, long_flye_medaka, long_unicycler, long_hybracter, short_unicycler, short_spades, hybrid_unicycler, hybrid_hybracter, hybrid_flye_medaka_pilon
	short_read1: Path?     // Illumina reads
	short_read2: Path?     // null for single-end reads
	long_reads: Path?      // ONT reads
	args: AssemblerArgs?
}

// An assembled genome, with the reads of its job (the input of GENOME_QC)
record AssembledGenome {
	meta: Record           // the meta of its job
	path: Path             // the assembly FASTA
	dir: Path              // the assembler outputs folder
	short_read1: Path?
	short_read2: Path?
	long_reads: Path?
}


// The assembly jobs of a sample: one per strategy of opts.name, named by it (assembler_name). r holds the meta of the
// sample (record(id)) and its reads, null when missing: read1, read2 (Illumina) and reads (ONT), e.g. the join of
// ShortReads and LongReads channels with remainder: true
def assemblyJobs(r: Record, opts: Record?) -> Iterable<AssemblyJob> {    // opts: AssemblerParams (typed Record: a record literal is not cast to a record type at run time)
	def names = (opts?.name ?: 'long_flye_medaka').tokenize(',').collect({ n -> n.trim() })
	return names.collect({ n -> record(
		meta: (r.meta as Record) + record(assembler_name: n),
		assembler: n,
		short_read1: r.read1 as Path,
		short_read2: r.read2 as Path,
		long_reads: r.reads as Path,
		args: opts?.args
	) })
}


workflow READS_ASSEMBLE {
	take:
	jobs: Channel<AssemblyJob>    // the assemblies to do (several strategies allowed, in any combination)

	main:
	def flye = FLYE(jobs.filter({ r -> r.assembler == 'long_flye' && r.long_reads }).map({ r -> r + record(args: r.args?.flye) }))
	def flye_medaka = LONG_FLYE_MEDAKA(jobs.filter({ r -> r.assembler == 'long_flye_medaka' && r.long_reads }))
	def flye_medaka_pilon = HYBRID_FLYE_MEDAKA_PILON(jobs.filter({ r -> r.assembler == 'hybrid_flye_medaka_pilon' && r.long_reads && r.short_read1 }))
	def spades = SPADES(jobs.filter({ r -> r.assembler == 'short_spades' && r.short_read1 }).map({ r -> r + record(long_reads: null, args: r.args?.spades) }))
	def unicycler = UNICYCLER(
		jobs.filter({ r -> r.assembler == 'long_unicycler' && r.long_reads }).map({ r -> r + record(short_read1: null, short_read2: null) })
			.mix(jobs.filter({ r -> r.assembler == 'short_unicycler' && r.short_read1 }).map({ r -> r + record(long_reads: null) }))
			.mix(jobs.filter({ r -> r.assembler == 'hybrid_unicycler' && r.long_reads && r.short_read1 }))
			.map({ r -> r + record(args: r.args?.unicycler) })
	)
	def hybracter = HYBRACTER(
		jobs.filter({ r -> r.assembler == 'long_hybracter' && r.long_reads }).map({ r -> r + record(short_read1: null, short_read2: null) })
			.mix(jobs.filter({ r -> r.assembler == 'hybrid_hybracter' && r.long_reads && r.short_read1 }))
			.map({ r -> r + record(args: r.args?.hybracter) })
	)

	def assembly_ch: Channel<Assembly> = flye
		.mix(flye_medaka)
		.mix(flye_medaka_pilon)
		.mix(spades)
		.mix(unicycler)
		.mix(hybracter)

	emit:
	// the assembly FASTA and folder, with the reads and the meta of its job (a single emit: the workflow returns this channel)
	genomes: Channel<AssembledGenome> = assembly_ch
		.join(jobs.map({ j -> record(meta: j.meta, short_read1: j.short_read1, short_read2: j.short_read2, long_reads: j.long_reads) }), by: 'meta')
}
