
#' Identifier of genes in a set of genomes
#'
#' @param genome Genomes of the genes.
#' @param id prokka ids of the genes (possibly the same in two genomes).
#' @return `<genome>:<id>`.
gene_uid <- function(genome, id) {
	str_c(genome, ":", id)
}

#' Core genes of genomes: clusters of reciprocal BLAST hits of the genes, with one gene in at least n_core genomes
#'
#' The genes are those of the reciprocal hits, each identified by `<genome>:<id>` (see gene_uid()) and linked to the
#' genes of its reciprocal hits. The genomes are those of the hits (a genome is paired with itself).
#'
#' @param hits A tibble of the reciprocal BLAST hits of the genes of each pair of genomes (itself included): the
#'   genomes of the query (query_genome) and of the subject (subject_genome), and the file of their hits (file), as
#'   written by the BLAST_RECIPROCAL_HITS subworkflow in the format of read_blast6().
#' @param n_core Minimum number of genomes sharing a core gene.
#' @return A list of two tibbles: `genes` (gene_id, genome, ID, length of the gene, cluster_id, is_core) and
#'   `clusters` (cluster_id, cluster_size, distinct_genome, mean_length, is_core), the clusters with more than one gene
#'   per genome being removed.
core_genes <- function(hits, n_core) {
	library(tidygraph)
	hits <- hits |>
		mutate(blast = map(file, read_blast6)) |>
		unnest(blast) |>
		mutate(query_gene = gene_uid(query_genome, query_id), subject_gene = gene_uid(subject_genome, subject_id))
	genes <- bind_rows(
		select(hits, gene_id = query_gene, genome = query_genome, ID = query_id, length = query_len),
		select(hits, gene_id = subject_gene, genome = subject_genome, ID = subject_id, length = subject_len)
	) |>
		distinct(gene_id, .keep_all = TRUE)

	# Clusters of genes linked by their reciprocal hits, keeping one paralog per genome in each cluster
	g <- tbl_graph(
		nodes = genes,
		edges = select(hits, from = query_gene, to = subject_gene, bit_score, pct_ident),
		directed = FALSE,
		node_key = "gene_id"
	) |>
		mutate(cluster_id = group_components()) |>
		mutate(score = centrality_closeness(weights = 100.1 - pct_ident)) |>
		group_by(cluster_id, genome) |>
		mutate(rk = row_number(score)) |>
		ungroup() |>
		filter(rk <= 1)

	genes <- as_tibble(g)
	clusters <- genes |>
		group_by(cluster_id) |>
		summarize(cluster_size = n_distinct(gene_id), distinct_genome = n_distinct(genome), mean_length = mean(length)) |>
		filter(cluster_size == distinct_genome) |>
		mutate(is_core = distinct_genome >= n_core)
	genes <- genes |>
		select(gene_id, genome, ID, length, cluster_id) |>
		left_join(select(clusters, cluster_id, is_core), by = "cluster_id") |>
		mutate(is_core = replace_na(is_core, FALSE))
	list(genes = genes, clusters = clusters)
}

#' Write the core genes, with their prokka product, and their sequences for pyMLST
#'
#' Writes:
#' - `core_genes.tsv`: the genes, with the prokka product of their protein;
#' - `pymlst_genomes/<genome>.fasta`: the core genes of each genome, for pyMLST;
#' - `pymlst_template.fasta`: the core genes of the genome with the most core genes, template of the pyMLST database.
#'
#' @param genes The genes, as returned by core_genes() (gene_id, genome, ..., is_core).
#' @param ffn A tibble of the genomes (genome) and of the nucleotide sequences of their genes (file: prokka.ffn, header:
#'   id and product).
core_genes_fasta <- function(genes, ffn) {
	stopifnot(!anyDuplicated(ffn$genome))
	seqs <- map2(ffn$genome, ffn$file, function(genome, file) {
		x <- Biostrings::readDNAStringSet(file)
		names(x) <- gene_uid(genome, names(x))    # <genome>:<id> <product>
		x
	}) |>
		do.call(what = c)
	products <- tibble(gene_id = str_extract(names(seqs), "^\\S+"), product = str_remove(names(seqs), "^\\S+\\s*"))
	names(seqs) <- products$gene_id
	genes |>
		left_join(products, by = "gene_id") |>
		relocate(product, .after = ID) |>
		write_tsv("core_genes.tsv")

	# Sequences of the core genes of each genome
	core <- filter(genes, is_core)
	dna <- split(seqs[core$gene_id], core$genome)
	dir.create("pymlst_genomes")
	for (genome in names(dna)) writeXStringSet(dna[[genome]], file.path("pymlst_genomes", paste0(genome, ".fasta")))
	template <- names(dna)[which.max(lengths(dna))]
	writeXStringSet(if (length(template)) dna[[template]] else DNAStringSet(), "pymlst_template.fasta")
}
