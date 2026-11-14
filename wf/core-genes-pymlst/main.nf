nextflow.enable.types = true

include { CORE_GENES_PYMLST } from './workflows/core-genes-pymlst/main.nf'
include { INPUT_FASTA; IdPath } from './subworkflows/input/main.nf'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'


// FASTA files of the genomes (one genome per file), and the options of CORE_GENES_PYMLST
params {
	fasta: String
	n_core: Integer                               // minimum number of genomes sharing a core gene (required)
	prokka_args: String = '--kingdom Bacteria'
	blastn_args: String = '-subject_besthit'      // alignment of the genes
	pymlst_args: String = '--coverage 0.85 --identity 0.85'    // `wgMLST add` options
}


workflow {
	main:
	// Validate parameters and print summary of supplied ones
	validateParameters()
	log.info(paramsSummaryLog(workflow))

	// the genomes: one per FASTA file, named by its sample id (gzipped files are decompressed)
	def genomes: Channel<IdPath> = INPUT_FASTA(params.fasta)
	def core = CORE_GENES_PYMLST(genomes, params.n_core, params.prokka_args, params.blastn_args, params.pymlst_args)

	publish:
	prokka = core.prokka
	core_genes = core.core_genes
	pymlst = core.pymlst
	report = core.report
}


output {
	prokka: Channel<Record> {
		path { r -> r.path >> "core-genes-pymlst/prokka/${r.meta.id}" }    // prokka of each genome
	}
	core_genes: Value<Record> {
		path { r ->
			r.core_genes    >> 'core-genes-pymlst/core_genes.tsv'
			r.core_graph    >> 'core-genes-pymlst/core_graph.rds'
			r.core_schema	>> 'core-genes-pymlst/core_schema.fasta'
		}
	}
	pymlst: Value<Record> {
		path { r ->
			r.dir >> 'core-genes-pymlst/pymlst'
		}
	}
	report: Value<Record> {
		path { r -> r.html >> 'core-genes-pymlst/report.html' }    // distance and allele heatmaps, minimum spanning tree
	}
}
