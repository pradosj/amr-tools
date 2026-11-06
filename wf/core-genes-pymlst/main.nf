nextflow.enable.types = true

include { CORE_GENES_PYMLST; Genome } from './workflows/core-genes-pymlst/main.nf'
include { INPUT_FASTA } from './subworkflows/input/main.nf'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'


// FASTA files of the genomes (one genome per file), and the options of CORE_GENES_PYMLST
params {
	fasta: String
	n_core: Integer                               // minimum number of genomes sharing a core gene (required)
	prokka_args: String = '--kingdom Bacteria'
	blastn_args: String = '-subject_besthit'      // alignment of the genes
	pymlst_args: String = '--coverage 0.85 --identity 0.85'    // `wgMLST add` options
}


// The genome named by each index of CORE_GENES_PYMLST (e.g. the pyMLST strains): a TSV with columns index, genome
process INDEX_TSV {
	input:
	text: String
	output:
	file('core_genes_index.tsv')
	exec:    // native task: the file is written by Nextflow, without a shell script
	task.workDir.resolve('core_genes_index.tsv').text = text
}


workflow {
	main:
	// Validate parameters and print summary of supplied ones
	validateParameters()
	log.info(paramsSummaryLog(workflow))

	// the genomes: one per FASTA file, named by its sample_id (gzipped files are decompressed)
	def genomes: Channel<Genome> = INPUT_FASTA(params.fasta).map({ r -> record(meta: r.meta, fasta: r.fasta) })
	def core = CORE_GENES_PYMLST(genomes, params.n_core, params.prokka_args, params.blastn_args, params.pymlst_args)

	// the genome of each index naming the genomes in the core genes and pyMLST results
	def index_tsv = INDEX_TSV(core.index.collect().map({ rs ->
		"index\tgenome\n" + rs.toList().toSorted({ r -> r.index }).collect({ r -> "${r.index}\t${r.meta.sample_id}\n" }).join('')
	}))

	publish:
	prokka = core.prokka
	core_genes = core.core_genes
	index = index_tsv
	pymlst = core.pymlst
}


output {
	prokka: Channel<Record> {
		path { r -> r.dir >> "prokka/${r.meta.sample_id}" }    // prokka of each genome
	}
	core_genes: Value<Record> {
		path { r ->
			r.genes_tsv >> 'core_genes.tsv'
			r.clusters_tsv >> 'core_clusters.tsv'
			r.template >> 'pymlst/template.fasta'
		}
	}
	index: Value<Path> {
		path { f -> f >> 'core_genes_index.tsv' }    // the genome of each index (pyMLST strain, genome column of the tables)
	}
	pymlst: Value<Record> {
		path { r ->
			r.db >> 'pymlst/pymlst.db'
			r.dir >> 'pymlst/results'
		}
	}
}
