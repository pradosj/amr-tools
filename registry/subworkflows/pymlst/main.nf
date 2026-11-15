nextflow.enable.types = true

include { PYMLST_WGMLST } from './modules/pymlst/wgmlst'
include { QUARTO_RENDER } from './modules/quarto/render'
include { MetaPath; IdPath } from './subworkflows/input/main.nf'


// A set of genomes to type on a schema with pyMLST
record PymlstInput {
	meta: Record?              // passed through to the outputs, unique to each set (e.g. the name of the set)
	schema: MetaPath           // FASTA of the genes of the schema
	genomes: List<IdPath>      // genes (or sequences) of each genome, the strains named by the id of their meta
	args: String?              // `wgMLST add` options (default: --coverage 0.85 --identity 0.85)
}


// pyMLST (wgMLST) typing of sets of genomes on a schema, and the HTML report of each set (assets/pymlst_report.qmd):
// distance and allele heatmaps, minimum spanning tree of the genomes (the genomes at distance 0 merged)
workflow PYMLST {
	take:
	inputs: Channel<PymlstInput>

	main:
	def typing_ch = PYMLST_WGMLST(inputs.map({ r -> r + record(args: r.args ?: '--coverage 0.85 --identity 0.85') }))
	def report_ch = QUARTO_RENDER(typing_ch.map({ r -> record(
		meta: r.meta,
		template: file("${moduleDir}/assets/pymlst_report.qmd"),
		extra: null,
		files: [r.dir],    // pymlst/ (distance.tsv, mlst.tsv)
		args: "-P 'input_dir:pymlst'"    // the files are staged next to the template
	) }))

	emit:
	pymlst = typing_ch.join(report_ch, by: 'meta')    // records with fields meta, db, dir (pyMLST outputs), html (report)
}
