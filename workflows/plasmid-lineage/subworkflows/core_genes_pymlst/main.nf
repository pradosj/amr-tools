nextflow.enable.types = true

include { PROKKA        } from './modules/tseemann/prokka'
include { BLASTP        } from './modules/blast/blastp'
include { PYMLST_WGMLST } from './modules/pymlst/wgmlst'


// The plasmids of a set, in one FASTA file
record PlasmidsFasta {
	meta: Map
	fasta: Path
}

// prokka annotation of the plasmids of a set
record ProkkaAnnotation {
	meta: Map
	dir: Path
}

// Core genes of a set of plasmids (see the CORE_GENES process)
record CoreGenes {
	meta: Map
	genes_tsv: Path
	clusters_tsv: Path
	template: Path
	genomes: List<Path>
	strains_tsv: Path
}

// pyMLST typing of the core genes of a set of plasmids
record PymlstTyping {
	meta: Map
	db: Path
	dir: Path    // distance.tsv, mlst.tsv, msa.fasta
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
	libs: List<Path>    // lib_plasmid_lineage.R, lib_core_gene.R

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
	source("lib/lib_plasmid_lineage.R")
	source("lib/lib_core_gene.R")
	core_genes("prokka", "prokka.self.blastp", n_core = ${n_core == null ? 'NA' : "${n_core}"})
	"""

	stub:
	"""
	mkdir pymlst_genomes
	touch core_genes.tsv core_clusters.tsv pymlst_template.fasta pymlst_genomes/plasmid.fasta pymlst_strains.tsv
	"""
}



// Core genes of a set of plasmids, typed with pyMLST: the prokka genes are clustered by reciprocal best BLASTP
// hits, and the clusters with one gene in at least n_core plasmids are the core genes
workflow CORE_GENES_PYMLST {
	take:
	plasmids: Channel<PlasmidsFasta>
	n_core: Integer?         // minimum number of plasmids sharing a core gene (null: 80% of the plasmids)
	prokka_args: String?
	blastp_args: String?     // alignment of the proteins
	pymlst_args: String?     // `wgMLST add` options

	main:
	def prokka_ch = PROKKA(plasmids.map({ r -> record(meta: r.meta, fasta: r.fasta, args: prokka_args) }))
	def blastp_ch = BLASTP(prokka_ch.map({ r -> record(
		meta: r.meta,
		query: r.dir.resolve('prokka.faa'),
		subject: null,
		args: blastp_args ?: '-subject_besthit'
	) }))
	def core_genes_ch = CORE_GENES(
		prokka_ch.join(blastp_ch, by: 'meta').map({ r -> record(meta: r.meta, prokka: r.dir, blast: r.tsv, n_core: n_core) }),
		[file("${moduleDir}/assets/lib_plasmid_lineage.R"), file("${moduleDir}/assets/lib_core_gene.R")]
	)
	def pymlst_ch = PYMLST_WGMLST(core_genes_ch.map({ r -> record(
		meta: r.meta,
		template: r.template,
		genomes: r.genomes,
		args: pymlst_args ?: '--coverage 0.85 --identity 0.85'
	) }))

	emit:
	prokka: Channel<ProkkaAnnotation> = prokka_ch
	core_genes: Channel<CoreGenes>    = core_genes_ch
	pymlst: Channel<PymlstTyping>     = pymlst_ch
}
