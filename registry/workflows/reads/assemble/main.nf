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

// Options of the assemblies (params.assembler of the pipelines, turned into AssemblyJob records by them: one job per strategy)
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
	asm_info: Path         // asm.json (see ASM_INFO): the strategy, the args of the assembly tools and the topology of the contigs
	short_read1: Path?
	short_read2: Path?
	long_reads: Path?
}


def assemblyTopology(strategy: String, dir: Path) -> List<Map> {
	if (strategy.contains('flye')) {
		def info = ['assembly_info.txt', '01_flye/assembly_info.txt', '01_flye_medaka/01_flye/assembly_info.txt']
			.collect({ f -> dir.resolve(f) })
			.find({ p -> p.exists() })
		if (!info) return []
		def rows = info.readLines().findAll({ l -> l }).collect({ l -> l.tokenize('\t') }).toList()
		return rows.subList(1, rows.size()).collect({ c -> [contig_id: c[0], topology: c[3] == 'Y' ? 'circular' : 'linear'] }).toList()
	}
	if (strategy.contains('unicycler')) {
		def fasta = dir.resolve('assembly.fasta')
		if (!fasta.exists()) return []
		return fasta.readLines().findAll({ l -> l.startsWith('>') }).collect({ l -> l.substring(1).tokenize(' ') })
			.collect({ h -> [contig_id: h[0], topology: h.contains('circular=true') ? 'circular' : 'linear'] }).toList()
	}
	if (strategy.contains('hybracter')) {
		def stats = dir.resolve('FINAL_OUTPUT/complete/sample_per_contig_stats.tsv')
		if (!stats.exists()) return []
		def rows = stats.readLines().findAll({ l -> l }).collect({ l -> l.tokenize('\t') }).toList()
		def c = rows[0].indexOf('circular')
		if (c < 0)
			return []
		return rows.subList(1, rows.size()).collect({ r -> [contig_id: r[0], topology: r[c] == 'True' ? 'circular' : 'linear'] }).toList()
	}
	return []
}


process ASM_INFO {
	input:
	record(meta: Record?, dir: Path, strategy: String, args: Record?)
	output:
	record(meta: meta, asm_info: file('asm.json'))
	exec:
	def info = [strategy: strategy, args: args, topology: assemblyTopology(strategy, dir)]
	task.workDir.resolve('asm.json').text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(info)) + '\n'
}


workflow READS_ASSEMBLE {
	take:
	jobs: Channel<AssemblyJob>

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

	def asm_info_ch = ASM_INFO(
		assembly_ch.join(jobs.map({ j -> record(meta: j.meta, strategy: j.assembler, args: j.args) }), by: 'meta')
	)

	emit:
	genomes: Channel<AssembledGenome> = assembly_ch
		.join(asm_info_ch, by: 'meta')
		.join(jobs.map({ j -> record(meta: j.meta, short_read1: j.short_read1, short_read2: j.short_read2, long_reads: j.long_reads) }), by: 'meta')
}
