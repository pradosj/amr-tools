nextflow.enable.types = true

include { PYMLST_WGMLST } from './modules/pymlst/wgmlst'
include { QUARTO_RENDER } from './modules/quarto/render'
include { MetaPath; IdPath } from './subworkflows/input/main.nf'


workflow PYMLST {
	take:
	schema: Value<MetaPath>
	genomes: Value<List<IdPath>>
	args: String?

	main:
	def typing_ch = PYMLST_WGMLST(channel.value(record(args: args ?: '--coverage 0.85 --identity 0.85')).combine(
		schema: schema,
		genomes: genomes
	))
	def report_ch = QUARTO_RENDER(typing_ch.map({ r -> record(
		meta: r.meta,
		template: file("${moduleDir}/assets/pymlst_report.qmd"),
		extra: null,
		files: [r.dir],
		args: "-P 'input_dir:pymlst'"
	) }))

	emit:
	pymlst = typing_ch    // records with fields meta (null), db, dir
	report = report_ch    // records with fields meta (null), html
}
