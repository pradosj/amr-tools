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

// Core genes from the reciprocal BLASTN hits of the genes of genomes and their sequences for pyMLST 
// Each gene is identified by the index of its genome
process CORE_GENES {
	container 'registry.gitlab.unige.ch/amr-genomics/rscript:v3'
	label 'rscript'
	memory 8.GB
	cpus 1
	time 30.m

	input:
	hits: List<GenomePairHits>
	genomes: List<GenomeGenes>
	n_core: Integer
	assets: Path    // the assets folder of the workflow (lib_io.R, lib_core_gene.R)

	stage:
	stageAs(hits.collect({ h -> h.tsv }), 'reciprocal_hits/?/*')
	stageAs(genomes.collect({ g -> g.ffn }), 'genes/?/*')
	stageAs(assets, 'assets')

	output:
	record(
	  core_graph: file('core_graph.rds'),
	  core_genes: file('core_genes.tsv'),
		core_template: file('core_template.fasta')
	)

	script:
	// the staged files with the index of their genomes, given to R as JSON
	def inputs = groovy.json.JsonOutput.toJson([
		hits: hits.collect({ h -> [query_genome: h.index_query, subject_genome: h.index_subject, file: "${h.tsv}"] }),
		ffn: genomes.collect({ g -> [genome: g.index, file: "${g.ffn}"] })
	])
	"""
	#!/usr/bin/env Rscript
	source("assets/lib_core_gene.R")
	inputs <- jsonlite::fromJSON('${inputs}')
	g <- core_genes(as_tibble(inputs[['hits']]),as_tibble(inputs[['ffn']]), n_core = ${n_core})
	saveRDS(g,'core_graph.rds')
	write_tsv(as_tibble(g,'nodes'),'core_genes.tsv')
	
	as_tibble(g,'nodes') %>% 
	  filter(is_core_template) %>% 
	  select(genome,gene_id,gene_uid) %>% 
	  left_join(inputs[['ffn']],by='genome') %>% 
	  group_by(genome,file) %>% 
	  group_map(~{read_fasta(.y\$file)[.x\$gene_id] |> setNames(.x\$gene_uid)}) %>% 
	  as("DNAStringSetList") %>% 
	  unlist() %>% 
	  Biostrings::writeXStringSet('core_template.fasta')
	"""

	stub:
	"""
	touch core_genes.tsv core_graph.rds core_template.fasta
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
	def blastn_ch = BLASTN(
	  prokka_ch.combine(db_ch)
	  .filter({ q, s -> q.meta != s.meta})
  	.map({ q, s -> record(
  		query: record(meta: q.meta, fasta: q.ffn),
  		subject: record(meta: s.meta, fasta: s.db),
  		args: blastn_args ?: '-subject_besthit'
  	)})
  )
	
  // pairing of reciprocal BLASTs and filtering to keep for each query the worst of bests reciprocal alignments
	def hits = BLAST_RECIPROCAL_HITS(blastn_ch, 0.75 as Float)
	
	// The genes of each genome with its index, sorted by index (the same order in each run for the cache)
	def indexed_genomes_collection = prokka_ch
		.join(index_ch, by: 'meta')
		.collect()
		.map({ rs -> rs.toList().toSorted({ r -> r.index }) })

	// Core genes: clusters of the reciprocal hits of each pair of genomes, and extract template sequences
	def core_genes_ch = CORE_GENES(
		hits
			.join(index_ch.map({ r -> record(meta_query: r.meta, index_query: r.index) }), by: 'meta_query')
			.join(index_ch.map({ r -> record(meta_subject: r.meta, index_subject: r.index) }), by: 'meta_subject')
			.collect()
			.map({ rs -> rs.toList().toSorted({ a, b -> a.index_query <=> b.index_query ?: a.index_subject <=> b.index_subject }) }), // the same order in each run for the cache
		indexed_genomes_collection,
		n_core,
		file("${moduleDir}/assets")
	)

	// pyMLST typing of the core genes: the strains are the indexes of the genomes (their position in genes)
	def pymlst_ch = PYMLST_WGMLST(core_genes_ch.combine(genomes: indexed_genomes_collection).map({ r -> record(
		meta: record(name: 'core_genes'),
		template: r.core_template,
		genomes: r.genomes.collect({ g -> g.ffn }).toList(),
		args: pymlst_args ?: '--coverage 0.85 --identity 0.85'
	) }))

	emit:
	index = index_ch              // records with fields meta (of the genome), index (naming the genome in core_genes and pymlst)
	prokka = prokka_ch            // records with fields meta (of the genome), dir, faa, ffn, gff, fna
	core_genes = core_genes_ch    // records with fields meta (of the genome), core_graph, core_genes, core_template
	pymlst = pymlst_ch            // records with fields meta (record(name: 'core_genes')), db, dir
}
