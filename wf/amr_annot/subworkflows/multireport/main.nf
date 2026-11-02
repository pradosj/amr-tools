nextflow.enable.types = true

include { RMD_RENDER        } from '../../modules/rmd/render'


// An entry of the AMR results folder: a file or folder, or a text written to a file
record NamedEntry {
	name: String      // path in the folder, may contain sub-folders
	path: Path?       // file or folder to put there...
	text: String?     // ...or text to write there
}


// Gather the AMR results into a single folder, each entry under its given name
process AMR_FOLDER {
	executor 'local'
	input:
	entries: List<NamedEntry>
	stage:
	stageAs(entries.findAll({ e -> e.path }).collect({ e -> e.path }), 'in/?/*')    // files may share the same name
	output:
	file('aggregated_folder')
	script:
	def cmds = entries.collect({ e ->
		def dest = "aggregated_folder/${e.name}"
		"mkdir -p \"\$(dirname '${dest}')\" && " + (e.path ? "mv '${e.path}' '${dest}'" : "echo '${e.text}' > '${dest}'")
	}).join('\n')
	"""
	mkdir aggregated_folder
	${cmds}
	"""
}


// Summary tables of the AMR results folder
process MULTITABLE {
	container 'registry.gitlab.unige.ch/amr-genomics/rscript:main'
	memory 8.GB
	cpus 2
	time 30.m
	input:
	db: Path     // AMR results folder
	lib: Path    // lib_typing.R
	output:
	record(
		xlsx: file('multitable.xlsx'),
		assemblies_tsv: file('multitable_assemblies.tsv'),
		contigs_tsv: file('multitable_contigs.tsv'),
		resistances_tsv: file('multitable_resistances.tsv'),
		plasmidfinder_hits_tsv: file('multitable_plasmidfinder_hits.tsv')
	)
	script:
	"""
	#!/usr/bin/env Rscript
	source("${lib}")
	db <- db_load("${db}")
	tbl <- list(
		assemblies = summarise_assembly(db),
		contigs = summarise_contigs(db),
		resistances = summarise_resistances(db),
		plasmidfinder_hits = summarise_plasmidfinder_hits(db)
	)
	tbl\$contigs <- left_join(tbl\$contigs,select(tbl\$assemblies,sample_id,assembler_name,species_name,mlst_type),by=c('sample_id','assembler_name'))
	openxlsx::write.xlsx(tbl,file="multitable.xlsx")
	readr::write_tsv(tbl\$assemblies,file="multitable_assemblies.tsv")
	readr::write_tsv(tbl\$contigs,file="multitable_contigs.tsv")
	readr::write_tsv(tbl\$resistances,file="multitable_resistances.tsv")
	readr::write_tsv(tbl\$plasmidfinder_hits,file="multitable_plasmidfinder_hits.tsv")
	"""
	stub:
	"""
	touch multitable.xlsx multitable_assemblies.tsv multitable_contigs.tsv multitable_resistances.tsv multitable_plasmidfinder_hits.tsv
	"""
}


// Gather the AMR results into a folder, then render the HTML report and the summary tables
workflow MULTIREPORT {
	take:
	entries: Channel<NamedEntry>

	main:
	def lib = file("${moduleDir}/assets/lib_typing.R")
	def folder_ch = AMR_FOLDER(
		entries.collect().flatMap({ es -> es.isEmpty() ? [] : [es.toList()] })  // no report without results
	)
	def html_ch = RMD_RENDER(folder_ch.map({ d -> record(
		template: file("${moduleDir}/assets/multireport.Rmd"),
		inputs: [d, lib],
		render_params: "indir='${d.name}'"
	) }))
	def table_ch = MULTITABLE(folder_ch, lib)

	emit:
	folder: Channel<Path> = folder_ch
	html: Channel<Path>   = html_ch
	tables: Channel<Record> = table_ch    // records with fields xlsx, assemblies_tsv, contigs_tsv, resistances_tsv, plasmidfinder_hits_tsv
}
