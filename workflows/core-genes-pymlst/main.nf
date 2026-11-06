nextflow.enable.types = true

include { PROKKA        } from './modules/tseemann/prokka'
include { MAKEBLASTDB   } from './modules/blast/makeblastdb'
include { BLASTP        } from './modules/blast/blastp'
include { PYMLST_WGMLST } from './modules/pymlst/wgmlst'


// The meta of a genome: its sample_id names it in the results (other fields allowed)
record GenomeMeta {
	sample_id: String
}

// A genome: its sequences in one FASTA file
record Genome {
	meta: GenomeMeta
	fasta: Path
}

// The prokka folder of a genome, named by the sample_id of the genome
record GenomeProkka {
	genome: String
	dir: Path
}

// The BLASTP of the proteins of a genome (query) against those of a genome (subject), named by their sample_id
record GenomePairBlast {
	query: String
	subject: String
	tsv: Path
}


// Core genes from the reciprocal BLASTP hits of the prokka proteins of genomes, and their sequences for pyMLST.
// Each gene is identified by the sample_id of its genome and its prokka id (prokka ids may be the same in two genomes).
process CORE_GENES {
	container 'registry.gitlab.unige.ch/amr-genomics/rscript:v3'
	label 'rscript'
	memory 8.GB
	cpus 1
	time 30.m

	input:
	record(meta: Map, prokka: List<GenomeProkka>, blast: List<GenomePairBlast>, n_core: Integer?)
	libs: List<Path>    // lib_io.R, lib_core_gene.R

	stage:
	stageAs(prokka.collect({ p -> p.dir }), 'prokka/?/*')
	stageAs(blast.collect({ b -> b.tsv }), 'blastp/?/*')
	stageAs(libs, 'lib/*')

	output:
	record(
		meta: meta,
		genes_tsv: file('core_genes.tsv'),          // genes with their genome and cluster, and whether the cluster is core
		clusters_tsv: file('core_clusters.tsv'),    // size and number of genomes of each cluster
		template: file('pymlst_template.fasta'),    // core genes of the genome with the most core genes
		genomes: files('pymlst_genomes/*.fasta').toSorted()    // core genes of each genome (pymlst_genomes/<sample_id>.fasta)
	)

	script:
	"""
	#!/usr/bin/env Rscript
	source("lib/lib_io.R")
	source("lib/lib_core_gene.R")
	prokka <- tibble(genome = ${rVector(prokka.collect({ p -> p.genome }))}, dir = ${rVector(prokka.collect({ p -> p.dir }))})
	blastp <- tibble(
		query_genome = ${rVector(blast.collect({ b -> b.query }))},
		subject_genome = ${rVector(blast.collect({ b -> b.subject }))},
		file = ${rVector(blast.collect({ b -> b.tsv }))}
	)
	core_genes(prokka, blastp, n_core = ${n_core == null ? 'NA' : "${n_core}"})
	"""

	stub:
	"""
	mkdir pymlst_genomes
	touch core_genes.tsv core_clusters.tsv pymlst_template.fasta pymlst_genomes/genome.fasta
	"""
}



// Core genes of genomes, typed with pyMLST: the prokka genes of the genomes are clustered by reciprocal BLASTP hits
// (BLASTP of each pair of genomes), and the clusters with one gene in at least n_core genomes are the core genes.
// The genomes are named by their sample_id (pyMLST strains).
workflow CORE_GENES_PYMLST {
	take:
	genomes: Channel<Genome>
	n_core: Integer?         // minimum number of genomes sharing a core gene (null: 80% of the genomes)
	prokka_args: String?
	blastp_args: String?     // alignment of the proteins
	pymlst_args: String?     // `wgMLST add` options

	main:
	// Annotate each genome, and make a BLAST database of its proteins
	def prokka_ch = PROKKA(genomes.map({ g -> record(meta: g.meta, fasta: g.fasta, args: prokka_args) }))
	def db_ch = MAKEBLASTDB(prokka_ch.map({ r -> record(meta: r.meta, fasta: r.faa, dbtype: 'prot', args: null) }))

	// BLASTP of the proteins of each genome (meta1) against those of each genome (meta2), itself included (paralogs)
	def blastp_ch = BLASTP(prokka_ch.combine(db_ch).map({ q, s -> record(
		meta: [meta1: q.meta, meta2: s.meta],
		query: q.faa,
		subject: null,
		db: s.db,
		args: blastp_args ?: '-subject_besthit'
	) }))

	// All the prokka folders and BLASTP results, named by the sample_id of their genomes
	def all_prokka_ch = prokka_ch
		.map({ r -> record(genome: (r.meta as GenomeMeta).sample_id, dir: r.dir) })
		.collect()
		.map({ rs -> rs.toList().toSorted({ r -> r.genome }) })
	def all_blastp_ch = blastp_ch
		.map({ r -> record(query: (r.meta.meta1 as GenomeMeta).sample_id, subject: (r.meta.meta2 as GenomeMeta).sample_id, tsv: r.tsv) })
		.collect()
		.map({ rs -> rs.toList().toSorted({ r -> "${r.query}\t${r.subject}" }) })

	def core_genes_ch = CORE_GENES(
		all_prokka_ch.combine(all_blastp_ch).map({ ps, bs -> record(meta: [:], prokka: ps, blast: bs, n_core: n_core) }),
		[file("${moduleDir}/assets/lib_io.R"), file("${moduleDir}/assets/lib_core_gene.R")]
	)
	def pymlst_ch = PYMLST_WGMLST(core_genes_ch.map({ r -> record(
		meta: r.meta,
		template: r.template,
		genomes: r.genomes,
		args: pymlst_args ?: '--coverage 0.85 --identity 0.85'
	) }))

	emit:
	prokka = prokka_ch            // records with fields meta (of the genome), dir, faa, gff, fna
	core_genes = core_genes_ch    // records with fields meta ([:]), genes_tsv, clusters_tsv, template, genomes
	pymlst = pymlst_ch            // records with fields meta ([:]), db, dir
}

// R character vector of values: c('x1', 'x2', ...)
def rVector(xs: Iterable) -> String {
	return "c(${xs.collect({ x -> "'${x}'" }).join(', ')})"
}
