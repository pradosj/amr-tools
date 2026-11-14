nextflow.enable.types = true

include { WG_BLASTN_DISTANCE } from './workflows/wg-blastn-distance/main.nf'
include { INPUT_FASTA; IdPath } from './subworkflows/input/main.nf'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'


// FASTA files of the genomes (one genome per file), and the options of WG_BLASTN_DISTANCE
params {
	fasta: String
	blastn_args: String = '-perc_identity 70 -evalue 1e-3 -subject_besthit -max_target_seqs 1000000'
	min_len: Integer = 1000
}


workflow {
	main:
	// Validate parameters and print summary of supplied ones
	validateParameters()
	log.info(paramsSummaryLog(workflow))

	// the genomes: one per FASTA file, named by its sample id (gzipped files are decompressed)
	def genomes: Channel<IdPath> = INPUT_FASTA(params.fasta)
	def distance_ch = WG_BLASTN_DISTANCE(genomes, params.blastn_args, params.min_len)

	publish:
	distance = distance_ch
}


output {
	distance: Value<Record> {
		path { r ->
			r.graph    >> 'wg-blastn-distance/distance_graph.rds'
		}
	}
}
