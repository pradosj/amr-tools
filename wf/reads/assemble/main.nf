nextflow.enable.types = true

include { READS_ASSEMBLE; AssemblerParams } from './workflows/assemble/main.nf'
include { INPUT_READS; ReadsParams                   } from './subworkflows/input/main.nf'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'


// Input reads given as glob patterns (the sample id is taken from the file names), and READS_ASSEMBLE options
params {
	reads: ReadsParams = record(ont: null, illumina: null)
	assembler: AssemblerParams = record(name: 'long_flye_medaka', args: record(flye: null, medaka: '--bacteria', spades: null, unicycler: null, hybracter: null))
	verbose: Boolean = false
}

workflow {
	main:
	// Validate parameters and print summary of supplied ones
	validateParameters()
	log.info(paramsSummaryLog(workflow))

	def reads = INPUT_READS(params.reads?.ont, params.reads?.illumina)
	// One assembly of each sample by each strategy of --assembler.name
	def assemblers = channel.fromList((params.assembler?.name ?: 'long_flye_medaka').tokenize(',').collect({ n -> n.trim() }))
	def assemblies = READS_ASSEMBLE(
		reads.short_reads
			.join(reads.long_reads.map({ r -> record(meta: r.meta, long_reads: r.path) }), by: 'meta', remainder: true)
			.combine(assemblers)
			.map({ r, a -> record(
				meta: r.meta + record(assembler_name: a),
				assembler: a,
				short_read1: r.read1,
				short_read2: r.read2,
				long_reads: r.long_reads,
				args: params.assembler?.args
			) })
	)

	publish:
	assemblies = assemblies
}

output {
	assemblies: Channel<Record> {
		mode "copy" //because link is problematic
		path { r ->
			r.path >> "assemblies/${r.meta.id}__${r.meta.assembler_name}.fasta"
			r.info >> "assemblies/${r.meta.id}__${r.meta.assembler_name}.asm.json"    // what the assembler tells of the contigs (topology), when it does
			r.dir >> (params.verbose ? "assemblies/${r.meta.id}__${r.meta.assembler_name}.output" : null)    // assembler folder, with --verbose
		}
	}
}
