nextflow.enable.types = true

include { PROKKA        } from './modules/tseemann/prokka'
include { MAKEBLASTDB   } from './modules/blast/makeblastdb'
include { BLASTN        } from './modules/blast/blastn'
include { PYMLST_WGMLST } from './modules/pymlst/wgmlst'
include { BLAST_RECIPROCAL_HITS } from './subworkflows/blast_reciprocal_hits/main.nf'

// A genome: its sequences in one FASTA file
record Genome {
	meta: Record
	fasta: Path
}

// The nucleotide sequences of the genes of a genome (prokka.ffn, header: id and product), named by the index of the
// genome (e.g. a PROKKA output joined with its index; other fields allowed)
record GenomeGenes {
	index: Integer
	ffn: Path
}

// The reciprocal hits of a pair of genomes (blast6 format, see BLAST_RECIPROCAL_HITS), named by the index of the
// genomes of their query and subject (e.g. a BLAST_RECIPROCAL_HITS output joined with the indexes; other fields allowed)
record GenomePairHits {
	index_query: Integer
	index_subject: Integer
	tsv: Path
}

// Core genes from the reciprocal BLASTN hits of the genes of genomes (see core_genes()), and their sequences for
// pyMLST (see core_genes_fasta()). Each gene is identified by the index of its genome and its id (ids may be the same
// in two genomes).
process CORE_GENES {
	container 'registry.gitlab.unige.ch/amr-genomics/rscript:v3'
	label 'rscript'
	memory 8.GB
	cpus 1
	time 30.m

	input:
	hits: List<GenomePairHits>
	ffn: List<GenomeGenes>
	n_core: Integer
	assets: Path    // the assets folder of the workflow (lib_io.R, lib_core_gene.R)

	stage:
	stageAs(hits.collect({ h -> h.tsv }), 'reciprocal_hits/?/*')
	stageAs(ffn.collect({ g -> g.ffn }), 'genes/?/*')
	stageAs(assets, 'assets')

	output:
	record(
		genes_tsv: file('core_genes.tsv'),          // genes with their genome, product and cluster, and whether the cluster is core
		clusters_tsv: file('core_clusters.tsv'),    // size and number of genomes of each cluster
		template: file('pymlst_template.fasta'),    // core genes of the genome with the most core genes
		genomes: files('pymlst_genomes/*.fasta').toSorted()    // core genes of each genome (pymlst_genomes/<index>.fasta)
	)

	script:
	// the staged files with the index of their genomes, given to R as JSON
	def inputs = groovy.json.JsonOutput.toJson([
		hits: hits.collect({ h -> [query_genome: h.index_query, subject_genome: h.index_subject, file: "${h.tsv}"] }),
		ffn: ffn.collect({ g -> [genome: g.index, file: "${g.ffn}"] })
	])
	"""
	#!/usr/bin/env Rscript
	source("assets/lib_io.R")
	source("assets/lib_core_gene.R")
	inputs <- jsonlite::fromJSON('${inputs}')
	core <- core_genes(as_tibble(inputs\$hits), n_core = ${n_core})
	write_tsv(core\$clusters, "core_clusters.tsv")
	core_genes_fasta(core\$genes, as_tibble(inputs\$ffn))
	"""

	stub:
	"""
	mkdir pymlst_genomes
	touch core_genes.tsv core_clusters.tsv pymlst_template.fasta pymlst_genomes/genome.fasta
	"""
}


// Core genes of genomes, typed with pyMLST: the prokka genes of the genomes are clustered by reciprocal BLASTN hits
// (BLASTN of the genes of each pair of genomes), and the clusters with one gene in at least n_core genomes are the core genes.
// The genomes are named by their index (from 1, in the order of the metas: see the index emit), e.g. the pyMLST strains.
workflow CORE_GENES_PYMLST {
	take:
	genomes: Channel<Genome>
	n_core: Integer          // minimum number of genomes sharing a core gene
	prokka_args: String?
	blastn_args: String?     // alignment of the genes
	pymlst_args: String?     // `wgMLST add` options

	main:
	// Index of each genome, naming it in the core genes and pyMLST results
	def index_ch = genomes
		.collect()
		.map({ gs -> gs.toList().collect({ g -> g.meta }).toSorted({ m -> "${m}" }) })
		.flatMap({ ms -> ms.withIndex().collect({ m, i -> record(meta: m, index: i + 1) }) })

	// Annotate each genome, and make a BLAST database of its genes
	def prokka_ch = PROKKA(genomes.map({ g -> record(meta: g.meta, fasta: g.fasta, args: prokka_args) }))
	def db_ch = MAKEBLASTDB(prokka_ch.map({ r -> record(meta: r.meta, fasta: r.ffn, dbtype: 'nucl', args: null) }))

	// BLASTN of the genes of each genome (query) against those of each genome (subject), itself included (paralogs)
	def blastn_ch = BLASTN(prokka_ch.combine(db_ch).map({ q, s -> record(
		query: record(meta: q.meta, path: q.ffn),
		subject: record(meta: s.meta, path: s.db),
		args: blastn_args ?: '-subject_besthit'
	) }))

	def hits = BLAST_RECIPROCAL_HITS(blastn_ch, 0.75 as Float)
	
	// Core genes: clusters of the reciprocal hits of each pair of genomes, and their sequences
	def core_genes_ch = CORE_GENES(
		hits
			.join(index_ch.map({ r -> record(meta_query: r.meta, index_query: r.index) }), by: 'meta_query')
			.join(index_ch.map({ r -> record(meta_subject: r.meta, index_subject: r.index) }), by: 'meta_subject')
			.collect()
			// the same order in each run for the cache
			.map({ rs -> rs.toList().toSorted({ a, b -> a.index_query <=> b.index_query ?: a.index_subject <=> b.index_subject }) }),
		prokka_ch
			.join(index_ch, by: 'meta')
			.collect()
			// the same order in each run for the cache
			.map({ rs -> rs.toList().toSorted({ r -> r.index }) }),
		n_core,
		file("${moduleDir}/assets")
	).map({ r -> r + record(meta: record(name: 'core_genes')) })

	// pyMLST typing of the core genes
	def pymlst_ch = PYMLST_WGMLST(core_genes_ch.map({ r -> record(
		meta: r.meta,
		template: r.template,
		genomes: r.genomes,
		args: pymlst_args ?: '--coverage 0.85 --identity 0.85'
	) }))

	emit:
	index = index_ch              // records with fields meta (of the genome), index (naming the genome in core_genes and pymlst)
	prokka = prokka_ch            // records with fields meta (of the genome), dir, faa, ffn, gff, fna
	core_genes = core_genes_ch    // records with fields meta (record(name: 'core_genes')), genes_tsv, clusters_tsv, template, genomes
	pymlst = pymlst_ch            // records with fields meta (record(name: 'core_genes')), db, dir
}
