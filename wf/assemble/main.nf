nextflow.enable.types = true

include { ASSEMBLE; AssemblerParams } from './workflows/assemble/main.nf'
include { INPUT_READS; ReadsParams                   } from './subworkflows/input/main.nf'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'


// Input reads given as glob patterns (the sample id is taken from the file names), and ASSEMBLE options
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
	def assemblies = ASSEMBLE(params.assembler, reads.short_reads, reads.long_reads.map({ r -> record(meta: r.meta, reads: r.path) }))

	publish:
	assemblies_fasta = assemblies
	assemblies_dir = assemblies
}

output {
	assemblies_fasta: Channel<Record> {
		mode "copy" //because link is problematic
		path { r -> r.fasta >> "assemblies/${r.meta.id}__${r.meta.assembler_name}.fasta" }
	}
	assemblies_dir: Channel<Record> {
		mode "copy" //because link is problematic
		enabled params.verbose
		path { r -> r.dir >> "assemblies/${r.meta.id}__${r.meta.assembler_name}.output" }
	}
}
