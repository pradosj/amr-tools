nextflow.enable.types = true

include { FLYE                          } from './modules/asm/flye'
include { LONG_FLYE_MEDAKA              } from './subworkflows/flye_medaka_pilon'
include { HYBRID_FLYE_MEDAKA_PILON      } from './subworkflows/flye_medaka_pilon'
include { UNICYCLER                     } from './modules/asm/unicycler'
include { SPADES                        } from './modules/asm/spades'
include { HYBRACTER                     } from './modules/asm/hybracter'
include { MetaPath                      } from './subworkflows/input/main.nf'


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


// An assembly to do: the reads of a sample, the strategy and the tools arguments
record AssemblyJob {
	meta: Record                // id (sample id) and assembler_name (names the assembly in the outputs, unique to each job)
	assembler: String           // assembly strategy: long_flye, long_flye_medaka, long_unicycler, long_hybracter, short_unicycler, short_spades, hybrid_unicycler, hybrid_hybracter, hybrid_flye_medaka_pilon
	short_reads: ShortReads?    // Illumina reads (their meta is not used), null when there are none
	long_reads: LongReads?      // ONT reads (their meta is not used), null when there are none
	args: AssemblerArgs?
}

// An assembled genome, with the reads it was assembled from
record AssembledGenome {
	meta: Record                // the meta of its job
	path: Path                  // the assembly FASTA
	short_reads: ShortReads?
	long_reads: LongReads?
}


// The assembly jobs of a sample: one per strategy of opts.name, named by it (assembler_name). r holds the meta of the
// sample (record(id)) and its reads, null when missing: read1, read2 (Illumina) and reads (ONT), e.g. the join of
// ShortReads and LongReads channels with remainder: true
def assemblyJobs(r: Record, opts: Record?) -> Iterable<AssemblyJob> {    // opts: AssemblerParams (typed Record: a record literal is not cast to a record type at run time)
	def short_reads: ShortReads? = null
	if (r.read1)
		short_reads = record(meta: null, read1: r.read1 as Path, read2: r.read2 as Path)
	def long_reads: LongReads? = null
	if (r.reads)
		long_reads = record(meta: null, reads: r.reads as Path)
	def names = (opts?.name ?: 'long_flye_medaka').tokenize(',').collect({ n -> n.trim() })
	return names.collect({ n -> record(
		meta: (r.meta as Record) + record(assembler_name: n),
		assembler: n,
		short_reads: short_reads,
		long_reads: long_reads,
		args: opts?.args
	) })
}


workflow READS_ASSEMBLE {
	take:
	jobs: Channel<AssemblyJob>    // the assemblies to do (several strategies allowed, in any combination)

	main:
	def flye = FLYE(jobs.filter({ r -> (r.assembler == 'long_flye') && (r.long_reads != null) }).map({ r -> record(meta: r.meta, long_reads: r.long_reads.reads, args: r.args?.flye) }))
	def flye_medaka = LONG_FLYE_MEDAKA(jobs.filter({ r -> (r.assembler == 'long_flye_medaka') && (r.long_reads != null) }).map({ r -> record(meta: r.meta, long_reads: r.long_reads.reads, args: r.args) }))
	def flye_medaka_pilon = HYBRID_FLYE_MEDAKA_PILON(jobs.filter({ r -> (r.assembler == 'hybrid_flye_medaka_pilon') && (r.long_reads != null) && (r.short_reads != null) }).map({ r -> record(meta: r.meta, short_read1: r.short_reads.read1, short_read2: r.short_reads.read2, long_reads: r.long_reads.reads, args: r.args) }))
	def spades = SPADES(jobs.filter({ r -> (r.assembler == 'short_spades') && (r.short_reads != null) }).map({ r -> record(meta: r.meta, short_read1: r.short_reads.read1, short_read2: r.short_reads.read2, args: r.args?.spades) }))
	def unicycler = UNICYCLER(
		jobs.filter({ r -> (r.assembler == 'long_unicycler') && (r.long_reads != null) }).map({ r -> record(meta: r.meta, short_read1: null, short_read2: null, long_reads: r.long_reads.reads, args: r.args?.unicycler) })
			.mix(jobs.filter({ r -> (r.assembler == 'short_unicycler') && (r.short_reads != null) }).map({ r -> record(meta: r.meta, short_read1: r.short_reads.read1, short_read2: r.short_reads.read2, long_reads: null, args: r.args?.unicycler) }))
			.mix(jobs.filter({ r -> (r.assembler == 'hybrid_unicycler') && (r.long_reads != null) && (r.short_reads != null) }).map({ r -> record(meta: r.meta, short_read1: r.short_reads.read1, short_read2: r.short_reads.read2, long_reads: r.long_reads.reads, args: r.args?.unicycler) }))
	)
	def hybracter = HYBRACTER(
		jobs.filter({ r -> (r.assembler == 'long_hybracter') && (r.long_reads != null) }).map({ r -> record(meta: r.meta, short_read1: null, short_read2: null, long_reads: r.long_reads.reads, args: r.args?.hybracter) })
			.mix(jobs.filter({ r -> (r.assembler == 'hybrid_hybracter') && (r.long_reads != null) && (r.short_reads != null) }).map({ r -> record(meta: r.meta, short_read1: r.short_reads.read1, short_read2: r.short_reads.read2, long_reads: r.long_reads.reads, args: r.args?.hybracter) }))
	)

	def fasta_ch: Channel<MetaPath> = flye.fasta
		.mix(flye_medaka.fasta)
		.mix(flye_medaka_pilon.fasta)
		.mix(spades.fasta)
		.mix(unicycler.fasta)
		.mix(hybracter.fasta)
	def outputs_ch: Channel<AssemblyOutputs> = flye.outputs
		.mix(flye_medaka.outputs)
		.mix(flye_medaka_pilon.outputs)
		.mix(spades.outputs)
		.mix(unicycler.outputs)
		.mix(hybracter.outputs)

	emit:
	// the assembly FASTA with the reads of its job, the meta of the job
	genomes: Channel<AssembledGenome> = fasta_ch
		.map({ r -> record(meta: r.meta, path: r.path) })
		.join(jobs.map({ j -> record(meta: j.meta, short_reads: j.short_reads, long_reads: j.long_reads) }), by: 'meta')
	outputs: Channel<AssemblyOutputs> = outputs_ch    // records with fields meta, fasta, dir (assembler outputs)
}
