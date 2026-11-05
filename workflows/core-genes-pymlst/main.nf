nextflow.enable.types = true

include { PROKKA        } from './modules/tseemann/prokka'
include { FASTA_CONCAT  } from './modules/fasta/concat'
include { BLASTP        } from './modules/blast/blastp'
include { PYMLST_WGMLST } from './modules/pymlst/wgmlst'


// A genome: its sequences in one FASTA file
record Genome {
	meta: Map
	fasta: Path
}


// Core genes from the reciprocal BLASTP hits of the prokka proteins of genomes, and their sequences for pyMLST
process CORE_GENES {
	container 'registry.gitlab.unige.ch/amr-genomics/rscript:v3'
	label 'rscript'
	memory 8.GB
	cpus 1
	time 30.m

	input:
	record(meta: Map, prokka: List<Path>, blast: Path, n_core: Integer?)    // prokka: the prokka folder of each genome
	libs: List<Path>    // lib_io.R, lib_core_gene.R

	stage:
	stageAs(prokka, 'prokka/?/*')
	stageAs(blast, 'prokka.self.blastp')
	stageAs(libs, 'lib/*')

	output:
	record(
		meta: meta,
		genes_tsv: file('core_genes.tsv'),          // genes with their genome and cluster, and whether the cluster is core
		clusters_tsv: file('core_clusters.tsv'),    // size and number of genomes of each cluster
		template: file('pymlst_template.fasta'),    // core genes of the genome with the most core genes
		genomes: files('pymlst_genomes/*.fasta').toSorted()    // core genes of each genome (pymlst_genomes/<genome>.fasta)
	)

	script:
	"""
	#!/usr/bin/env Rscript
	source("lib/lib_io.R")
	source("lib/lib_core_gene.R")
	core_genes(Sys.glob("prokka/*/*"), "prokka.self.blastp", n_core = ${n_core == null ? 'NA' : "${n_core}"})
	"""

	stub:
	"""
	mkdir pymlst_genomes
	touch core_genes.tsv core_clusters.tsv pymlst_template.fasta pymlst_genomes/g1.fasta
	"""
}



// Core genes of genomes, typed with pyMLST: the prokka genes of the genomes are clustered by reciprocal best BLASTP
// hits, and the clusters with one gene in at least n_core genomes are the core genes.
// The i-th genome (from 1, in the order of their meta) is named g<i>: the prefix of its prokka genes (--locustag)
// and its pyMLST strain.
workflow CORE_GENES_PYMLST {
	take:
	genomes: Channel<Genome>
	n_core: Integer?         // minimum number of genomes sharing a core gene (null: 80% of the genomes)
	prokka_args: String?
	blastp_args: String?     // alignment of the proteins
	pymlst_args: String?     // `wgMLST add` options

	main:
	// Name each genome g<i>: meta1 is the meta of the genome, meta2 its index i and its locustag g<i>
	def named_ch = genomes
		.collect()
		.flatMap({ gs -> gs.toList().toSorted({ g -> "${g.meta}" }).withIndex().collect({
			g, i -> record(meta: [meta1: g.meta, meta2: [index: i + 1, locustag: 'g' + (i + 1)]], fasta: g.fasta)
		}) })

	// Annotate each genome, its genes named <locustag>_<n>
	def prokka_ch = PROKKA(named_ch.map({ 
			r -> record(meta: r.meta, fasta: r.fasta, args: "${prokka_args ?: '--kingdom Bacteria'} --locustag ${r.meta.meta2.locustag}") 
		}))
	def all_prokka_ch = prokka_ch
		.collect()
		.map({ rs -> rs.toList().toSorted({ r -> r.meta.meta2.index }) })

	// Reciprocal alignment of the proteins of all the genomes
	def faa_ch = FASTA_CONCAT(all_prokka_ch.map({ rs -> record(meta: [:], fasta: rs.collect({ r -> r.faa }).toList()) }))
	def blastp_ch = BLASTP(faa_ch.map({ r -> record(
		meta: r.meta,
		query: r.fasta,
		subject: null,
		args: blastp_args ?: '-subject_besthit'
	) }))
	def core_genes_ch = CORE_GENES(
		all_prokka_ch.combine(blastp_ch).map({ rs, b -> record(meta: [:], prokka: rs.collect({ r -> r.dir }).toList(), blast: b.tsv, n_core: n_core) }),
		[file("${moduleDir}/assets/lib_io.R"), file("${moduleDir}/assets/lib_core_gene.R")]
	)
	def pymlst_ch = PYMLST_WGMLST(core_genes_ch.map({ r -> record(
		meta: r.meta,
		template: r.template,
		genomes: r.genomes,
		args: pymlst_args ?: '--coverage 0.85 --identity 0.85'
	) }))

	emit:
	prokka = prokka_ch            // records with fields meta ([meta1: meta of the genome, meta2: [index: i, locustag: g<i>]]), dir, faa, gff, fna
	core_genes = core_genes_ch    // records with fields meta ([:]), genes_tsv, clusters_tsv, template, genomes
	pymlst = pymlst_ch            // records with fields meta ([:]), db, dir
}
