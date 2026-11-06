# Functions of the CORE_GENES_PYMLST workflow (sourced after lib_io.R, whose functions it uses)

suppressPackageStartupMessages(library(GenomicRanges))

#' Minimum number of plasmids to call a gene or a region core
#'
#' @param n_core The given minimum number of plasmids, or NA.
#' @param n_plasmids Number of plasmids.
#' @return n_core, or 80% of the plasmids when n_core is NA.
core_count <- function(n_core, n_plasmids) {
	if (is.null(n_core) || is.na(n_core)) ceiling(0.8 * n_plasmids) else as.integer(n_core)
}

#' Core genes of genomes: clusters of reciprocal BLASTP hits of the prokka proteins, with one gene in at least n_core
#' genomes
#'
#' Each gene is identified by `<genome>:<prokka id>`, prokka ids being possibly the same in two genomes.
#'
#' Writes:
#' - `core_genes.tsv`: the genes of the clusters (gene_id, genome, prokka ID, prokka gene name), with their cluster
#'   and whether the cluster is core;
#' - `core_clusters.tsv`: the size and number of genomes of each cluster;
#' - `pymlst_genomes/<genome>.fasta`: the core genes of each genome, for pyMLST;
#' - `pymlst_template.fasta`: the core genes of the genome with the most core genes, template of the pyMLST database.
#'
#' @param prokka A tibble of the genomes (genome) and of their prokka output folder (dir: prokka.gff, prokka.ffn).
#' @param blastp A tibble of the BLASTP of the prokka proteins of each genome (query_genome) against those of each genome
#'   (subject_genome, itself included): file, as read by read_blast6().
#' @param n_core Minimum number of genomes sharing a core gene (NA: see core_count()).
#' @param min_cover Minimum coverage of the query protein by its alignment.
core_genes <- function(prokka, blastp, n_core = NA, min_cover = 0.75) {
	library(tidygraph)
	stopifnot(!anyDuplicated(prokka$genome))
	gene_id <- function(genome, id) str_c(genome, ":", id)
	gff <- map2(prokka$genome, prokka$dir, function(genome, dir) {
		x <- rtracklayer::import.gff3(file.path(dir, "prokka.gff"))
		x <- x[x$type == "CDS"]
		gene <- if (is.null(x$gene)) rep(NA_character_, length(x)) else as.character(x$gene)
		tibble(gene_id = gene_id(genome, x$ID), genome = genome, ID = x$ID, gene = gene, width = width(x))
	}) |>
		bind_rows()
	n_core <- core_count(n_core, nrow(prokka))

	# Reciprocal hits, keeping the worst direction of each pair of genes
	z <- pmap(blastp, function(query_genome, subject_genome, file) {
		read_blast6(file) |>
			mutate(query_id = gene_id(query_genome, query_id), subject_id = gene_id(subject_genome, subject_id))
	}) |>
		bind_rows() |>
		filter(align_len / query_len > min_cover) |>
		group_by(query_id, subject_id) |>
		slice_max(bit_score, n = 1, with_ties = FALSE) |>
		ungroup() |>
		filter(str_c(query_id, subject_id, sep = "|") %in% str_c(subject_id, query_id, sep = "|")) |>
		filter(query_id != subject_id) |>
		mutate(undirected_edge_id = str_c(pmin(query_id, subject_id), pmax(query_id, subject_id), sep = "|")) |>
		group_by(undirected_edge_id) |>
		slice_min(bit_score, n = 1, with_ties = FALSE) |>
		ungroup()

	# Clusters of genes, keeping one paralog per genome in each cluster
	g <- tbl_graph(
		nodes = select(gff, gene_id, genome, ID, gene, width),
		edges = select(z, from = query_id, to = subject_id, bit_score, pct_ident),
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
		summarize(cluster_size = n_distinct(gene_id), distinct_genome = n_distinct(genome), mean_width = mean(width)) |>
		filter(cluster_size == distinct_genome) |>
		mutate(is_core = distinct_genome >= n_core)
	genes <- genes |>
		select(gene_id, genome, ID, gene, width, cluster_id) |>
		left_join(select(clusters, cluster_id, is_core), by = "cluster_id") |>
		mutate(is_core = replace_na(is_core, FALSE))
	write_tsv(genes, "core_genes.tsv")
	write_tsv(clusters, "core_clusters.tsv")

	# Sequences of the core genes of each genome
	ffn <- map2(prokka$genome, prokka$dir, function(genome, dir) {
		x <- read_fasta(file.path(dir, "prokka.ffn"))
		names(x) <- gene_id(genome, names(x))
		x
	}) |>
		do.call(what = c)
	core <- filter(genes, is_core)
	dna <- split(ffn[core$gene_id], core$genome)
	dir.create("pymlst_genomes")
	for (genome in names(dna)) writeXStringSet(dna[[genome]], file.path("pymlst_genomes", paste0(genome, ".fasta")))
	template <- names(dna)[which.max(lengths(dna))]
	writeXStringSet(if (length(template)) dna[[template]] else DNAStringSet(), "pymlst_template.fasta")
}
