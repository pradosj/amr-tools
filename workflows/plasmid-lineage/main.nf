nextflow.enable.types = true

include { BLASTN as BLASTN_SELF } from './modules/blast/blastn'
include { QUARTO_RENDER         } from './modules/quarto/render'
include { AMR_ANNOT; AmrParams  } from './workflows/amr-annot/main.nf'
include { CORE_GENES_PYMLST     } from './workflows/core-genes-pymlst/main.nf'


// Options of PLASMID_LINEAGE (params.lineage). Fields are nullable: defaults are used when null.
record PlasmidLineageParams {
	core_genome_count: Integer?    // minimum number of plasmids sharing a core gene or core region (required)
	blastn_args: String?           // alignment of the plasmids
	genes_blastn_args: String?     // alignment of the genes, for the core genes
	prokka_args: String?
	pymlst_args: String?           // `wgMLST add` options
}

// A plasmid: its FASTA file (a genome of CORE_GENES_PYMLST and of AMR_ANNOT)
record Genome {
	meta: Record
	fasta: Path
}


// Gather the plasmid sequences in one FASTA (adding the final newline missing in some files)
process CONCAT_FASTA {
	label 'local'

	input:
	record(meta: Record, fasta: List<Path>)

	stage:
	stageAs(fasta, 'in/?/*')    // files may share the same name

	output:
	record(meta: meta, fasta: file('plasmids.fasta'))

	script:
	"""
	awk 1 ${fasta.collect({ f -> "'${f}'" }).join(' ')} > plasmids.fasta
	"""
}



// The plasmid named by each index of CORE_GENES_PYMLST (e.g. the pyMLST strains): a TSV with columns index, genome
// (sample_id of the plasmid) and contig_id (the id of each of its sequences)
process INDEX_TSV {
	input:
	text: String
	output:
	file('core_genes_index.tsv')
	exec:    // native task: the file is written by Nextflow, without a shell script
	task.workDir.resolve('core_genes_index.tsv').text = text
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
	source("assets/lib_plasmid_lineage.R")
	find_align_seq("plasmids.fasta", "plasmids.self.blastn", "align_seq.fasta")
	"""

	stub:
	"""
	touch align_seq.fasta
	"""
}


// Core genome: regions of the plasmids shared by most plasmids
process CORE_GENOME {
	container 'registry.gitlab.unige.ch/amr-genomics/rscript:v3'
	label 'rscript'
	memory 8.GB
	cpus 1
	time 30.m

	input:
	record(meta: Record, fasta: Path, blast: Path, n_core: Integer)
	assets: Path

	stage:
	stageAs(fasta, 'plasmids.fasta')
	stageAs(blast, 'plasmids.self.blastn')
	stageAs(assets, 'assets')

	output:
	record(
		meta: meta,
		tsv: file('core_genome.tsv')    // regions with the number of plasmids sharing them
	)

	script:
	"""
	#!/usr/bin/env Rscript
	source("assets/lib_io.R")
	source("assets/lib_plasmid_lineage.R")
	core_genome("plasmids.fasta", "plasmids.self.blastn", n_core = ${n_core})
	"""

	stub:
	"""
	touch core_genome.tsv
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


// Lineage of a set of plasmids (each a genome: a FASTA file): align the plasmids, annotate them, find their core
// genes (typed with pyMLST) and their core genome, and render an HTML report
workflow PLASMID_LINEAGE {
	take:
	opts: PlasmidLineageParams
	amr_opts: AmrParams
	amrfinderplus_db: Value<Path?>    // database of AMR_ANNOT (null when AMRFinderPlus is skipped)
	name: String                      // name of the set of plasmids (meta of its results: record(sample_id: name))
	plasmids: Channel<Genome>         // the plasmids, named by the sample_id of their meta

	main:
	def assets = file("${moduleDir}/assets")
	if (opts.core_genome_count == null) {
		error("PLASMID_LINEAGE: the minimum number of plasmids sharing a core gene or region (--lineage.core_genome_count) is required")
	}
	def n_core = opts.core_genome_count as Integer
	def blastn_args = opts.blastn_args ?: '-perc_identity 70 -evalue 1e-3 -subject_besthit'
	def set_meta = record(sample_id: name)

	// Align the plasmids against each other (all their sequences in one FASTA)
	def fasta_ch = CONCAT_FASTA(plasmids
		.map({ r -> r.fasta })
		.collect()
		.flatMap({ fs -> fs.isEmpty() ? [] : [record(meta: set_meta, fasta: fs.toList().toSorted())] })
	)
	def self_ch = BLASTN_SELF(fasta_ch.map({ r -> record(query: record(meta: r.meta, path: r.fasta), subject: null, args: blastn_args) }))
		.map({ r -> record(meta: r.meta_query, tsv: r.tsv) })

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

	// Core genes: prokka genes with reciprocal BLASTN hits in most plasmids, typed with pyMLST. The plasmids are named
	// by their index in the results (core_genes_index.tsv: the plasmid and the sequences of each index)
	def core = CORE_GENES_PYMLST(plasmids, n_core, opts.prokka_args, opts.genes_blastn_args, opts.pymlst_args)
	def index_tsv = INDEX_TSV(core.index
		.join(plasmids, by: 'meta')
		.collect()
		.map({ rs -> "index\tgenome\tcontig_id\n" + rs.toList().toSorted({ r -> r.index }).collect({ r ->
			fastaIds(r.fasta).collect({ id -> "${r.index}\t${r.meta.sample_id}\t${id}\n" }).join('')
		}).join('') })
	)
	def core_genes_ch = core.core_genes.combine(index_tsv: index_tsv).map({ r -> record(
		meta: set_meta,
		genes_tsv: r.genes_tsv,
		clusters_tsv: r.clusters_tsv,
		index_tsv: r.index_tsv,
		template: r.template,
		genomes: r.genomes
	) })
	def pymlst_ch = core.pymlst.map({ r -> record(meta: set_meta, db: r.db, dir: r.dir) })

	// Core genome: regions shared by most plasmids
	def core_genome_ch = CORE_GENOME(
		fasta_ch.join(self_ch, by: 'meta').map({ r -> record(meta: r.meta, fasta: r.fasta, blast: r.tsv, n_core: n_core) }),
		assets
	)

	// HTML report of the set of plasmids, from its files gathered in a folder under explicit names
	def report_data_ch = REPORT_FOLDER(
		fasta_ch.map({ r -> record(meta: r.meta, plasmids: r.fasta) })
			.join(self_ch.map({ r -> record(meta: r.meta, self_blast: r.tsv) }), by: 'meta')
			.join(core_genome_ch.map({ r -> record(meta: r.meta, core_genome_tsv: r.tsv) }), by: 'meta')
			.combine(core_genes: core_genes_ch)
			.combine(pymlst: pymlst_ch)
			.combine(amr.multiqc)
			.map({ r, m -> record(meta: r.meta, entries: [
				record(name: 'plasmids.fasta', path: r.plasmids),
				record(name: 'plasmids.self.blastn', path: r.self_blast),
				record(name: 'core_genes.tsv', path: r.core_genes.genes_tsv),
				record(name: 'core_clusters.tsv', path: r.core_genes.clusters_tsv),
				record(name: 'core_genes_index.tsv', path: r.core_genes.index_tsv),
				record(name: 'pymlst', path: r.pymlst.dir),
				record(name: 'core_genome.tsv', path: r.core_genome_tsv),
				record(name: 'multitable_resistances.tsv', path: m.resistances_tsv as Path),
				record(name: 'multitable_plasmidfinder_hits.tsv', path: m.plasmidfinder_hits_tsv as Path)
			]) })
	)
	def report_ch = QUARTO_RENDER(report_data_ch.map({ r -> record(
		meta: r.meta,
		template: file("${moduleDir}/assets/plasmid_lineage.qmd"),
		extra: assets,
		files: [r.dir],
		args: "-P 'name:${r.meta.sample_id}' -P 'core_genome_count:${n_core}'"
	) }))

	emit:
	prokka: Channel<Record>          = core.prokka        // records with fields meta (of the plasmid), dir, faa, ffn, gff, fna
	core_genes: Value<Record>        = core_genes_ch      // records with fields meta (of the set), genes_tsv, clusters_tsv, index_tsv, template, genomes
	pymlst: Value<Record>            = pymlst_ch          // records with fields meta (of the set), db, dir
	core_genome: Channel<Record>     = core_genome_ch     // records with fields meta (of the set), tsv
	report: Channel<Record>          = report_ch          // records with fields meta (of the set), html
	amr_assemblies: Channel<Record>  = amr.assemblies_amr    // records with fields meta (of the plasmid), fasta, fai, dir, rds
	amr_multiqc: Channel<Record>     = amr.multiqc           // records with fields rds, html, xlsx, assemblies_tsv, contigs_tsv, resistances_tsv, plasmidfinder_hits_tsv
}


// Ids of the sequences of a FASTA file: first word of their header
def fastaIds(fasta: Path) -> List<String> {
	return fasta.text.readLines().findAll({ l -> l.startsWith('>') }).collect({ l -> l.substring(1).tokenize(' \t')[0] }).toList()
}
