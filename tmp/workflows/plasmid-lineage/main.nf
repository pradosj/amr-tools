nextflow.enable.types = true

include { FASTA_CAT             } from './modules/fasta/cat'
include { MAKEBLASTDB           } from './modules/aln/blast/makeblastdb'
include { BLASTN as BLASTN_SELF } from './modules/aln/blast/blastn'
include { QUARTO_RENDER         } from './modules/r/quarto/render'
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
	name: String                      // name of the set of plasmids (meta of its results: record(id: name))
	plasmids: Channel<Genome>         // the plasmids

	main:
	def assets = file("${moduleDir}/assets")
	def blastn_args = opts.blastn_args ?: '-perc_identity 70 -evalue 1e-3 -subject_besthit'
	def set_meta = record(id: name)

	// Align the plasmids against each other (all their sequences in one FASTA)
	def fasta_ch = FASTA_CAT(plasmids
		.map({ r -> record(meta: null, path: r.fasta) })    // no meta: the sequence names are kept
		.collect()
		.flatMap({ rs -> rs.isEmpty() ? [] : [record(meta: set_meta, genomes: rs.toList().toSorted({ r -> r.path }))] })
	).fasta.map({ r -> record(meta: r.meta, fasta: r.path) })
	def self_ch = BLASTN_SELF(
		fasta_ch
			.join(MAKEBLASTDB(fasta_ch.map({ r -> record(meta: r.meta, path: r.fasta) })).map({ d -> record(meta: d.meta, db: d) }), by: 'meta')
			.map({ r -> record(meta: r.meta, path: r.fasta, ref: r.db, args: blastn_args) })
	)

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
			.join(self_ch.map({ r -> record(meta: r.meta, self_blast: r.path) }), by: 'meta')
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
		args: "-P 'name:${r.meta.id}'"
	) }))

	emit:
	report: Channel<Record>          = report_ch          // records with fields meta (of the set), html
	amr_assemblies: Channel<Record>  = amr.assemblies_amr    // records with fields meta (of the plasmid), fasta, fai, dir, rds
	amr_multiqc: Channel<Record>     = amr.multiqc           // records with fields rds, html, xlsx, assemblies_tsv, contigs_tsv, resistances_tsv, plasmidfinder_hits_tsv
}
