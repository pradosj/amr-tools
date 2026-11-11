nextflow.enable.types = true

include { MAKEBLASTDB } from './modules/blast/makeblastdb'
include { BLASTN as BLASTN_SELF } from './modules/blast/blastn'
include { QUARTO_RENDER         } from './modules/quarto/render'
include { AMR_ANNOT; AmrParams  } from './workflows/amr-annot/main.nf'


// Options of PLASMID_LINEAGE (params.lineage). Fields are nullable: defaults are used when null.
record PlasmidLineageParams {
	blastn_args: String?           // alignment of the plasmids
}

// A plasmid: its FASTA file (a genome of AMR_ANNOT)
record Genome {
	meta: Record
	fasta: Path
}




// Longest region shared by most plasmids (empty when there is none): the reference used by the pipeline to rotate the plasmids
process FIND_ALIGN_SEQ {
	container 'registry.gitlab.unige.ch/amr-genomics/rscript:v3'
	label 'rscript'
	memory 8.GB
	cpus 1
	time 30.m

	input:
	record(meta: Record, fasta: Path, blast: Path)
	assets: Path

	stage:
	stageAs(fasta, 'plasmids.fasta')
	stageAs(blast, 'plasmids.self.blastn')
	stageAs(assets, 'assets')

	output:
	record(meta: meta, fasta: file('align_seq.fasta'))

	script:
	"""
	#!/usr/bin/env Rscript
	source("assets/lib_io.R")
	source("assets/lib_blast.R")
	source("assets/lib_plasmid_lineage.R")
	find_align_seq("plasmids.fasta", "plasmids.self.blastn", "align_seq.fasta")
	"""

	stub:
	"""
	touch align_seq.fasta
	"""
}


// A file or folder of the report data, and its name in the report_data folder
record NamedFile {
	name: String
	path: Path
}

// Gather the files of the report in the report_data folder, each under its given name
process REPORT_FOLDER {
	label 'local'

	input:
	record(meta: Record, entries: List<NamedFile>)

	stage:
	stageAs(entries.collect({ e -> e.path }), 'in/?/*')    // files may share the same name

	output:
	record(meta: meta, dir: file('report_data'))

	script:
	def cmds = entries.collect({ e -> "cp -P '${e.path}' 'report_data/${e.name}'" }).join('\n')    // cp: a staged file can be given under several names
	"""
	mkdir report_data
	${cmds}
	"""
}

			

// Lineage of a set of plasmids (each a genome: a FASTA file): align the plasmids, annotate them, and render an HTML
// report
workflow PLASMID_LINEAGE {
	take:
	opts: PlasmidLineageParams
	amr_opts: AmrParams
	amrfinderplus_db: Value<Path?>    // database of AMR_ANNOT (null when AMRFinderPlus is skipped)
	name: String                      // name of the set of plasmids (meta of its results: record(sample_id: name))
	plasmids: Channel<Genome>         // the plasmids

	main:
	def assets = file("${moduleDir}/assets")
	def blastn_args = opts.blastn_args ?: '-perc_identity 70 -evalue 1e-3 -subject_besthit'
	def set_meta = record(sample_id: name)

	// Align the plasmids against each other (all their sequences in one FASTA)
	def fasta_ch = CONCAT_FASTA(plasmids
		.map({ r -> r.fasta })
		.collect()
		.flatMap({ fs -> fs.isEmpty() ? [] : [record(meta: set_meta, fasta: fs.toList().toSorted())] })
	)
	def self_ch = BLASTN_SELF(
		fasta_ch
			.join(MAKEBLASTDB(fasta_ch.map({ r -> record(meta: r.meta, fasta: r.fasta, dbtype: 'nucl', args: null) })), by: 'meta')
			.map({ r -> record(query: record(meta: r.meta, fasta: r.fasta), subject_db: record(meta: r.meta, db: r.db), args: blastn_args) })
	)
		.map({ r -> record(meta: r.meta_query, blast6: r.blast6) })

	// AMR and plasmid typing annotations of each plasmid, without orgfinder (and so without cgemlst, whose scheme
	// depends on the organism found by orgfinder)
	def amr = AMR_ANNOT(
		amr_opts + record(orgfinder: record(skip: true, args: null), cgemlst: record(skip: true, args: null)),
		null,
		amrfinderplus_db,
		plasmids.map({ r -> record(meta: r.meta + record(assembler_name: 'none'), fasta: r.fasta, gff: null, info: null) }),
		channel.empty(),
		channel.empty()
	)

	// HTML report of the set of plasmids, from its files gathered in a folder under explicit names
	def report_data_ch = REPORT_FOLDER(
		fasta_ch.map({ r -> record(meta: r.meta, plasmids: r.fasta) })
			.join(self_ch.map({ r -> record(meta: r.meta, self_blast: r.blast6) }), by: 'meta')
			.combine(amr.multiqc)
			.map({ r, m -> record(meta: r.meta, entries: [
				record(name: 'plasmids.fasta', path: r.plasmids),
				record(name: 'plasmids.self.blastn', path: r.self_blast),
				record(name: 'multitable_resistances.tsv', path: m.resistances_tsv as Path),
				record(name: 'multitable_plasmidfinder_hits.tsv', path: m.plasmidfinder_hits_tsv as Path)
			]) })
	)
	def report_ch = QUARTO_RENDER(report_data_ch.map({ r -> record(
		meta: r.meta,
		template: file("${moduleDir}/assets/plasmid_lineage.qmd"),
		extra: assets,
		files: [r.dir],
		args: "-P 'name:${r.meta.sample_id}'"
	) }))

	emit:
	report: Channel<Record>          = report_ch          // records with fields meta (of the set), html
	amr_assemblies: Channel<Record>  = amr.assemblies_amr    // records with fields meta (of the plasmid), fasta, fai, dir, rds
	amr_multiqc: Channel<Record>     = amr.multiqc           // records with fields rds, html, xlsx, assemblies_tsv, contigs_tsv, resistances_tsv, plasmidfinder_hits_tsv
}
