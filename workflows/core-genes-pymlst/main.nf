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

// The index of a genome (from 1, in the order of the metas), naming it in the core genes and pyMLST results
record GenomeIndex {
	meta: Record
	index: Integer
}

// The nucleotide sequences of the genes of a genome (prokka.ffn, header: id and product), named by the index of the
// genome
record GenomeGenes {
	genome: Integer
	ffn: Path
}

// The reciprocal hits of a pair of genomes (blast6 format, see BLAST_RECIPROCAL_HITS), named by the index of the
// genomes of their query and subject
record GenomePairHits {
	query: Integer
	subject: Integer
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
	libs: List<Path>    // lib_io.R, lib_core_gene.R

	stage:
	stageAs(hits.collect({ h -> h.tsv }), 'reciprocal_hits/?/*')
	stageAs(ffn.collect({ g -> g.ffn }), 'genes/?/*')
	stageAs(libs, 'lib/*')

	output:
	record(
		genes_tsv: file('core_genes.tsv'),          // genes with their genome, product and cluster, and whether the cluster is core
		clusters_tsv: file('core_clusters.tsv'),    // size and number of genomes of each cluster
		template: file('pymlst_template.fasta'),    // core genes of the genome with the most core genes
		genomes: files('pymlst_genomes/*.fasta').toSorted()    // core genes of each genome (pymlst_genomes/<index>.fasta)
	)

	script:
	"""
	#!/usr/bin/env Rscript
	source("lib/lib_io.R")
	source("lib/lib_core_gene.R")
	hits <- tibble(
		query_genome = ${rVector(hits.collect({ h -> h.query }))},
		subject_genome = ${rVector(hits.collect({ h -> h.subject }))},
		file = ${rVector(hits.collect({ h -> h.tsv }))}
	)
	ffn <- tibble(genome = ${rVector(ffn.collect({ g -> g.genome }))}, file = ${rVector(ffn.collect({ g -> g.ffn }))})
	core <- core_genes(hits, n_core = ${n_core})
	write_tsv(core\$clusters, "core_clusters.tsv")
	core_genes_fasta(core\$genes, ffn)
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
	def metas_ch = genomes.collect().map({ gs -> gs.toList().collect({ g -> g.meta }).toSorted({ m -> "${m}" }) })
	def index_ch = metas_ch.flatMap({ ms -> ms.withIndex().collect({ m, i -> record(meta: m, index: i + 1) }) })

	// Annotate each genome, and make a BLAST database of its genes
	def prokka_ch = PROKKA(genomes.map({ g -> record(meta: g.meta, fasta: g.fasta, args: prokka_args) }))
	def db_ch = MAKEBLASTDB(prokka_ch.map({ r -> record(meta: r.meta, fasta: r.ffn, dbtype: 'nucl', args: null) }))

	// BLASTN of the genes of each genome (query) against those of each genome (subject), itself included (paralogs)
	def blastn_ch = BLASTN(prokka_ch.combine(db_ch).map({ q, s -> record(
		query: record(meta: q.meta, path: q.ffn),
		subject: record(meta: s.meta, path: s.db),
		args: blastn_args ?: '-subject_besthit'
	) }))

	// Core genes: clusters of the reciprocal hits of each pair of genomes, and their sequences
	def hits = BLAST_RECIPROCAL_HITS(blastn_ch, 0.75 as Float)
	def all_hits_ch = hits
		.collect()
		.map({ rs -> record(hits: rs.toList()) })
		.combine(metas: metas_ch)
		.map({ r -> r.hits
			.collect({ h -> record(query: r.metas.indexOf(h.meta_query) + 1, subject: r.metas.indexOf(h.meta_subject) + 1, tsv: h.tsv) })
			.toSorted({ h -> h.query * 1000000 + h.subject })
		})
	def all_ffn_ch = prokka_ch
		.collect()
		.map({ rs -> record(prokka: rs.toList()) })
		.combine(metas: metas_ch)
		.map({ r -> r.prokka
			.collect({ p -> record(genome: r.metas.indexOf(p.meta) + 1, ffn: p.ffn) })
			.toSorted({ g -> g.genome })
		})
	def core_genes_ch = CORE_GENES(all_hits_ch, all_ffn_ch, n_core, [file("${moduleDir}/assets/lib_io.R"), file("${moduleDir}/assets/lib_core_gene.R")])
		.map({ r -> r + record(meta: record(name: 'core_genes')) })

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

// R character vector of values: c('x1', 'x2', ...)
def rVector(xs: Iterable) -> String {
	return "c(${xs.collect({ x -> "'${x}'" }).join(', ')})"
}
