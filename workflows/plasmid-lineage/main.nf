nextflow.enable.types = true

include { BLASTN as BLASTN_SELF; BLASTN as BLASTN_ROTATED } from './modules/blast/blastn'
include { BLASTP                } from './modules/blast/blastp'
include { PROKKA                } from './modules/tseemann/prokka'
include { PYMLST_WGMLST         } from './modules/pymlst/wgmlst'
include { FASTANI               } from './modules/fastani'
include { QUARTO_RENDER         } from './modules/quarto/render'
include { AMR_ANNOT; AmrParams  } from './workflows/amr-annot/main.nf'
include { ROTATE_GENOMES        } from './workflows/rotate-genomes/main.nf'


// Options of PLASMID_LINEAGE (params.lineage). Fields are nullable: defaults are used when null.
record PlasmidLineageParams {
	core_genome_count: Integer?    // minimum number of plasmids sharing a core gene or core region (default: 80% of the plasmids)
	blastn_args: String?           // alignment of the plasmids
	blastp_args: String?           // alignment of the proteins, for the core genes
	prokka_args: String?
	pymlst_args: String?           // `wgMLST add` options
	fastani_args: String?          // ANI of the core genomes
}

// A set of plasmids to compare: FASTA files with one or several plasmid sequences
record PlasmidSet {
	meta: Map           // [sample_id: name of the set]
	fasta: List<Path>
}


// Gather the plasmid sequences in one FASTA (adding the final newline missing in some files)
process CONCAT_FASTA {
	label 'local'

	input:
	record(meta: Map, fasta: List<Path>)

	stage:
	stageAs(fasta, 'in/?/*')    // files may share the same name

	output:
	record(meta: meta, fasta: file('plasmids.fasta'))

	script:
	"""
	awk 1 in/*/* > plasmids.fasta
	"""
}


// Longest region shared by most plasmids (empty when there is none): the origin of the rotated plasmids
process FIND_ALIGN_SEQ {
	container 'registry.gitlab.unige.ch/amr-genomics/rscript:v3'
	label 'rscript'
	memory 8.GB
	cpus 1
	time 30.m

	input:
	record(meta: Map, fasta: Path, blast: Path)
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
	source("assets/lib_plasmid_lineage.R")
	find_align_seq("plasmids.fasta", "plasmids.self.blastn", "align_seq.fasta")
	"""

	stub:
	"""
	touch align_seq.fasta
	"""
}


// Core genes from the reciprocal BLASTP hits of the prokka proteins, and their sequences for pyMLST
process CORE_GENES {
	container 'registry.gitlab.unige.ch/amr-genomics/rscript:v3'
	label 'rscript'
	memory 8.GB
	cpus 1
	time 30.m

	input:
	record(meta: Map, prokka: Path, blast: Path, n_core: Integer?)
	assets: Path

	stage:
	stageAs(prokka, 'prokka')
	stageAs(blast, 'prokka.self.blastp')
	stageAs(assets, 'assets')

	output:
	record(
		meta: meta,
		genes_tsv: file('core_genes.tsv'),          // genes with their cluster, and whether the cluster is core
		clusters_tsv: file('core_clusters.tsv'),    // size and number of plasmids of each cluster
		template: file('pymlst_template.fasta'),    // core genes of the plasmid with the most core genes
		genomes: files('pymlst_genomes/*.fasta').toSorted(),    // core genes of each plasmid
		strains_tsv: file('pymlst_strains.tsv')     // plasmid of each pyMLST strain
	)

	script:
	"""
	#!/usr/bin/env Rscript
	source("assets/lib_plasmid_lineage.R")
	core_genes("prokka", "prokka.self.blastp", n_core = ${n_core == null ? 'NA' : "${n_core}"})
	"""

	stub:
	"""
	mkdir pymlst_genomes
	touch core_genes.tsv core_clusters.tsv pymlst_template.fasta pymlst_genomes/plasmid.fasta pymlst_strains.tsv
	"""
}


// Core genome: regions of the rotated plasmids shared by most plasmids, and the regions of each plasmid for FastANI
process CORE_GENOME {
	container 'registry.gitlab.unige.ch/amr-genomics/rscript:v3'
	label 'rscript'
	memory 8.GB
	cpus 1
	time 30.m

	input:
	record(meta: Map, fasta: Path, blast: Path, n_core: Integer?)
	assets: Path

	stage:
	stageAs(fasta, 'rotated.fasta')
	stageAs(blast, 'rotated.self.blastn')
	stageAs(assets, 'assets')

	output:
	record(
		meta: meta,
		tsv: file('core_genome.tsv'),                            // regions with the number of plasmids sharing them
		genomes: files('fastani_genomes/*.fasta').toSorted(),    // regions of each plasmid
		genomes_tsv: file('fastani_genomes.tsv')                 // plasmid of each FastANI genome
	)

	script:
	"""
	#!/usr/bin/env Rscript
	source("assets/lib_plasmid_lineage.R")
	core_genome("rotated.fasta", "rotated.self.blastn", n_core = ${n_core == null ? 'NA' : "${n_core}"})
	"""

	stub:
	"""
	mkdir fastani_genomes
	touch core_genome.tsv fastani_genomes/plasmid.fasta fastani_genomes.tsv
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
	record(meta: Map, entries: List<NamedFile>)

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


// Lineage of a set of plasmids: rotate the plasmids on a shared region, annotate them, find their core
// genes (typed with pyMLST) and their core genome (compared with FastANI), and render an HTML report
workflow PLASMID_LINEAGE {
	take:
	opts: PlasmidLineageParams
	amr_opts: AmrParams
	org_db: Value<Path?>              // databases of AMR_ANNOT (null when the tool is skipped)
	amrfinderplus_db: Value<Path?>
	plasmids: Channel<PlasmidSet>

	main:
	def assets = file("${moduleDir}/assets")
	def n_core = opts.core_genome_count
	def blastn_args = opts.blastn_args ?: '-perc_identity 70 -evalue 1e-3 -subject_besthit'

	// Align the plasmids against each other
	def fasta_ch = CONCAT_FASTA(plasmids)
	def self_ch = BLASTN_SELF(fasta_ch.map({ r -> record(meta: r.meta, query: r.fasta, subject: null, args: blastn_args) }))

	// Rotate the plasmids to start on the longest region shared by most of them, and align them again
	def align_seq_ch = FIND_ALIGN_SEQ(
		fasta_ch.join(self_ch, by: 'meta').map({ r -> record(meta: r.meta, fasta: r.fasta, blast: r.tsv) }),
		assets
	)
	def rotated_ch = ROTATE_GENOMES(
		fasta_ch,
		align_seq_ch.map({ r -> record(meta: r.meta, fasta: r.fasta, blast: null) }),
		blastn_args
	).rotated
	def rotated_self_ch = BLASTN_ROTATED(rotated_ch.map({ r -> record(meta: r.meta, query: r.fasta, subject: null, args: blastn_args) }))

	// AMR and plasmid typing annotations of the rotated plasmids (one multiqc for all the sets)
	def amr = AMR_ANNOT(
		amr_opts,
		org_db,
		amrfinderplus_db,
		rotated_ch.map({ r -> record(meta: r.meta + [assembler_name: 'none'], fasta: r.fasta, gff: null, info: null) }),
		channel.empty(),
		channel.empty()
	)

	// Core genes: prokka genes with reciprocal BLASTP hits in most plasmids, typed with pyMLST
	def prokka_ch = PROKKA(rotated_ch.map({ r -> r + record(args: opts.prokka_args) }))
	def blastp_ch = BLASTP(prokka_ch.map({ r -> record(
		meta: r.meta,
		query: (r.dir as Path).resolve('prokka.faa'),
		subject: null,
		args: opts.blastp_args ?: '-subject_besthit'
	) }))
	def core_genes_ch = CORE_GENES(
		prokka_ch.join(blastp_ch, by: 'meta').map({ r -> record(meta: r.meta, prokka: r.dir, blast: r.tsv, n_core: n_core) }),
		assets
	)
	def pymlst_ch = PYMLST_WGMLST(core_genes_ch.map({ r -> record(
		meta: r.meta,
		template: r.template,
		genomes: r.genomes,
		args: opts.pymlst_args ?: '--coverage 0.85 --identity 0.85'
	) }))

	// Core genome: regions shared by most plasmids, compared with FastANI
	def core_genome_ch = CORE_GENOME(
		rotated_ch.join(rotated_self_ch, by: 'meta').map({ r -> record(meta: r.meta, fasta: r.fasta, blast: r.tsv, n_core: n_core) }),
		assets
	)
	def fastani_ch = FASTANI(core_genome_ch.map({ r -> record(
		meta: r.meta,
		query: r.genomes,
		refs: [],
		args: opts.fastani_args ?: '--fragLen 500'
	) }))

	// HTML report of each set of plasmids, from its files gathered in a folder under explicit names
	def report_data_ch = REPORT_FOLDER(
		fasta_ch.map({ r -> record(meta: r.meta, plasmids: r.fasta) })
			.join(self_ch.map({ r -> record(meta: r.meta, self_blast: r.tsv) }), by: 'meta')
			.join(rotated_ch.map({ r -> record(meta: r.meta, rotated: r.fasta) }), by: 'meta')
			.join(rotated_self_ch.map({ r -> record(meta: r.meta, rotated_blast: r.tsv) }), by: 'meta')
			.join(prokka_ch, by: 'meta')
			.join(core_genes_ch.map({ r -> record(meta: r.meta, genes_tsv: r.genes_tsv, clusters_tsv: r.clusters_tsv, strains_tsv: r.strains_tsv) }), by: 'meta')
			.join(pymlst_ch.map({ r -> record(meta: r.meta, pymlst: r.dir) }), by: 'meta')
			.join(core_genome_ch.map({ r -> record(meta: r.meta, core_genome_tsv: r.tsv, genomes_tsv: r.genomes_tsv) }), by: 'meta')
			.join(fastani_ch.map({ r -> record(meta: r.meta, fastani: r.tsv) }), by: 'meta')
			.combine(amr.multiqc)
			.map({ r, m -> record(meta: r.meta, entries: [
				record(name: 'plasmids.fasta', path: r.plasmids),
				record(name: 'plasmids.self.blastn', path: r.self_blast),
				record(name: 'rotated.fasta', path: r.rotated),
				record(name: 'rotated.self.blastn', path: r.rotated_blast),
				record(name: 'prokka', path: r.dir),
				record(name: 'core_genes.tsv', path: r.genes_tsv),
				record(name: 'core_clusters.tsv', path: r.clusters_tsv),
				record(name: 'pymlst_strains.tsv', path: r.strains_tsv),
				record(name: 'pymlst', path: r.pymlst),
				record(name: 'core_genome.tsv', path: r.core_genome_tsv),
				record(name: 'fastani_genomes.tsv', path: r.genomes_tsv),
				record(name: 'fastani.tsv', path: r.fastani),
				record(name: 'multitable_resistances.tsv', path: m.resistances_tsv as Path),
				record(name: 'multitable_plasmidfinder_hits.tsv', path: m.plasmidfinder_hits_tsv as Path)
			]) })
	)
	def report_ch = QUARTO_RENDER(report_data_ch.map({ r -> record(
		meta: r.meta,
		template: file("${moduleDir}/assets/plasmid_lineage.qmd"),
		extra: assets,
		files: [r.dir],
		args: "-P 'name:${r.meta.sample_id}'" + (n_core != null ? " -P 'core_genome_count:${n_core}'" : '')
	) }))

	emit:
	rotated: Channel<Record>         = rotated_ch         // records with fields meta, fasta
	prokka: Channel<Record>          = prokka_ch          // records with fields meta, dir
	core_genes: Channel<Record>      = core_genes_ch      // records with fields meta, genes_tsv, clusters_tsv, template, genomes, strains_tsv
	pymlst: Channel<Record>          = pymlst_ch          // records with fields meta, db, dir
	core_genome: Channel<Record>     = core_genome_ch     // records with fields meta, tsv, genomes, genomes_tsv
	fastani: Channel<Record>         = fastani_ch         // records with fields meta, tsv
	report: Channel<Record>          = report_ch          // records with fields meta, html
	amr_assemblies: Channel<Record>  = amr.assemblies_amr    // records with fields meta, fasta, fai, dir, rds
	amr_multiqc: Channel<Record>     = amr.multiqc           // records with fields rds, html, xlsx, assemblies_tsv, contigs_tsv, resistances_tsv, plasmidfinder_hits_tsv
}
