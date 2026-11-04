nextflow.enable.types = true

include { PROKKA        } from './modules/tseemann/prokka'
include { BLASTP        } from './modules/blast/blastp'
include { PYMLST_WGMLST } from './modules/pymlst/wgmlst'


// A genome: a set of sequences in one FASTA file (e.g. the plasmids of a set), whose core genes are the genes found
// in most of its sequences
record Genome {
	meta: Map
	fasta: Path
	gff: Path?
	info: Path?
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
	libs: List<Path>    // lib_io.R, lib_core_gene.R

	stage:
	stageAs(prokka, 'prokka')
	stageAs(blast, 'prokka.self.blastp')
	stageAs(libs, 'lib/*')

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
	source("lib/lib_io.R")
	source("lib/lib_core_gene.R")
	core_genes("prokka", "prokka.self.blastp", n_core = ${n_core == null ? 'NA' : "${n_core}"})
	"""

	stub:
	"""
	mkdir pymlst_genomes
	touch core_genes.tsv core_clusters.tsv pymlst_template.fasta pymlst_genomes/plasmid.fasta pymlst_strains.tsv
	"""
}



// Core genes of genomes, typed with pyMLST: in each genome (a set of sequences, e.g. plasmids), the prokka genes are
// clustered by reciprocal best BLASTP hits, and the clusters with one gene in at least n_core sequences are the core genes
workflow CORE_GENES_PYMLST {
	take:
	genomes: Channel<Genome>    // each genome is a set of sequences (e.g. plasmids): its core genes are shared by most of them
	n_core: Integer?         // minimum number of plasmids sharing a core gene (null: 80% of the plasmids)
	prokka_args: String?
	blastp_args: String?     // alignment of the proteins
	pymlst_args: String?     // `wgMLST add` options

	main:
	def prokka_ch = PROKKA(genomes.map({ r -> record(meta: r.meta, fasta: r.fasta, args: prokka_args) }))
	def blastp_ch = BLASTP(prokka_ch.map({ r -> record(
		meta: r.meta,
		query: r.faa,
		args: blastp_args ?: '-subject_besthit'
	) }))
	def core_genes_ch = CORE_GENES(
		prokka_ch.join(blastp_ch, by: 'meta').map({ r -> record(meta: r.meta, prokka: r.dir, blast: r.tsv, n_core: n_core) }),
		[file("${moduleDir}/assets/lib_io.R"), file("${moduleDir}/assets/lib_core_gene.R")]
	)
	def pymlst_ch = PYMLST_WGMLST(core_genes_ch.map({ r -> record(
		meta: r.meta,
		template: r.template,
		genomes: r.genomes,
		args: pymlst_args ?: '--coverage 0.85 --identity 0.85'
	) }))

	emit:
	prokka = prokka_ch            // records with fields meta, dir, faa, gff, fna
	core_genes = core_genes_ch    // records with fields meta, genes_tsv, clusters_tsv, template, genomes, strains_tsv
	pymlst = pymlst_ch            // records with fields meta, db, dir
}
