nextflow.enable.types = true

include { PROKKA        } from './modules/tseemann/prokka'
include { MAKEBLASTDB   } from './modules/blast/makeblastdb'
include { BLASTP        } from './modules/blast/blastp'
include { PYMLST_WGMLST } from './modules/pymlst/wgmlst'
include { BLAST_RECIPROCAL_HITS } from './subworkflows/blast_reciprocal_hits/main.nf'


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

// The reciprocal hits of a pair of genomes (blast6 format, see BLAST_RECIPROCAL_HITS), named by the sample_id of the
// genomes of their query and subject
record GenomePairHits {
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
	record(meta: Record, prokka: List<GenomeProkka>, hits: List<GenomePairHits>, n_core: Integer)
	libs: List<Path>    // lib_io.R, lib_core_gene.R

	stage:
	stageAs(prokka.collect({ p -> p.dir }), 'prokka/?/*')
	stageAs(hits.collect({ h -> h.tsv }), 'reciprocal_hits/?/*')
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
	hits <- tibble(
		query_genome = ${rVector(hits.collect({ h -> h.query }))},
		subject_genome = ${rVector(hits.collect({ h -> h.subject }))},
		file = ${rVector(hits.collect({ h -> h.tsv }))}
	)
	core_genes(hits, n_core = ${n_core}) |>
		write_core_genes(prokka)
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
	n_core: Integer          // minimum number of genomes sharing a core gene
	prokka_args: String?
	blastp_args: String?     // alignment of the proteins
	pymlst_args: String?     // `wgMLST add` options

	main:
	// Annotate each genome, and make a BLAST database of its proteins
	def prokka_ch = PROKKA(genomes.map({ g -> record(meta: g.meta, fasta: g.fasta, args: prokka_args) }))
	def db_ch = MAKEBLASTDB(prokka_ch.map({ r -> record(meta: r.meta, fasta: r.faa, dbtype: 'prot', args: null) }))

	// BLASTP of the proteins of each genome (query) against those of each genome (subject), itself included (paralogs)
	def blastp_ch = BLASTP(prokka_ch.combine(db_ch).map({ q, s -> record(
		query: record(meta: q.meta, path: q.faa),
		subject: record(meta: s.meta, path: s.db),
		args: blastp_args ?: '-subject_besthit'
	) }))

	// Reciprocal hits of each pair of genomes
	def hits = BLAST_RECIPROCAL_HITS(blastp_ch, null)

	// All the prokka folders, named by the sample_id of their genomes
	def all_prokka_ch = prokka_ch
		.map({ r -> record(genome: r.meta.sample_id, dir: r.dir) })
		.collect()
		.map({ rs -> rs.toList().toSorted({ r -> r.genome }) })

	// The reciprocal hits of all the pairs of genomes
	def all_hits_ch = hits
		.map({ r -> record(query: r.meta_query.sample_id as String, subject: r.meta_subject.sample_id as String, tsv: r.tsv) })
		.collect()
		.map({ rs -> rs.toList().toSorted({ r -> "${r.query}\t${r.subject}" }) })

	def core_genes_ch = CORE_GENES(
		all_prokka_ch
			.map({ ps -> record(prokka: ps) })
			.combine(hits: all_hits_ch)
			.map({ r -> record(meta: record(name: 'core_genes'), prokka: r.prokka, hits: r.hits, n_core: n_core) }),
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
	core_genes = core_genes_ch    // records with fields meta (record(name: 'core_genes')), genes_tsv, clusters_tsv, template, genomes
	pymlst = pymlst_ch            // records with fields meta (record(name: 'core_genes')), db, dir
}

// R character vector of values: c('x1', 'x2', ...)
def rVector(xs: Iterable) -> String {
	return "c(${xs.collect({ x -> "'${x}'" }).join(', ')})"
}
