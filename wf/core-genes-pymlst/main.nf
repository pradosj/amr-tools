nextflow.enable.types = true

include { CORE_GENES_PYMLST; Genome } from './workflows/core-genes-pymlst/main.nf'
include { INPUT_FASTA } from './subworkflows/input/main.nf'
include { QUARTO_RENDER } from './modules/quarto/render'
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

	// the genomes: one per FASTA file, named by its sample_id (gzipped files are decompressed)
	def genomes: Channel<Genome> = INPUT_FASTA(params.fasta).map({ r -> record(meta: record(id: r.meta.sample_id), fasta: r.fasta) })
	def core = CORE_GENES_PYMLST(genomes, params.n_core, params.prokka_args, params.blastn_args, params.pymlst_args)

	// HTML report: distance and allele heatmaps, minimum spanning tree of the pyMLST results
	def report_ch = QUARTO_RENDER(core.pymlst.combine(genes: core.core_genes.map({ r -> r.core_genes })).map({ r -> record(
		meta: r.meta,
		template: file("${projectDir}/assets/core_genes_pymlst.qmd"),
		extra: null,
		files: [r.dir, r.genes],    // pymlst/ (distance.tsv, mlst.tsv), core_genes.tsv
		args: "-P 'n_core:${params.n_core}' -P 'input_dir:./'"    // the files are staged next to the template
	) }))

	publish:
	prokka = core.prokka
	core_genes = core.core_genes
	pymlst = core.pymlst
	report = report_ch
}


output {
	prokka: Channel<Record> {
		path { r -> r.dir >> "core-genes-pymlst/prokka/${r.meta.id}" }    // prokka of each genome
	}
	core_genes: Value<Record> {
		path { r ->
			r.core_genes    >> 'core-genes-pymlst/core_genes.tsv'
			r.core_graph    >> 'core-genes-pymlst/core_graph.rds'
			r.core_template	>> 'core-genes-pymlst/core_template.fasta'
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
