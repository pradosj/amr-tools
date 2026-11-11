source("assets/lib_io.R")
source("assets/lib_blast.R")


#' Core genes of genomes: clusters of reciprocal BLAST hits of the genes, with one gene in at least n_core genomes
#'
#' The genes are those of the FASTA file of the genes of all the genomes, each named `<genome>:<id>` (gene_uid), as
#' written by FASTA_CAT with the id of the genome (which may contain ":", the gene id does not). Each gene is linked to the genes of other genomes of its
#' reciprocal hits (the hits between genes of the same genome are ignored).
#'
#' @param blast Self BLAST (BLASTN) of the genes FASTA file, in the format of read_blast6().
#' @param genes_fasta FASTA file of the genes of all the genomes, named `<genome>:<id>`.
#' @param n_core Minimum number of genomes sharing a core gene.
#' @param min_cover Minimum coverage of the query gene by its alignment (see blast_reciprocal_hits()).
#' @return A tbl_graph: the genes as nodes (gene_uid, genome, gene_id, cluster_id, score, is_core_template), a single
#'   gene per genome and cluster (the most central), linked by their reciprocal hits.
core_genes <- function(blast, genes_fasta, n_core, min_cover = 0.75) {
	genes <- tibble(gene_uid = names(read_fasta(genes_fasta))) |>
		mutate(genome = str_remove(gene_uid, ":[^:]*$"), gene_id = str_extract(gene_uid, "[^:]*$"))    # the genome id may contain ":"

	# Reciprocal hits between genes of distinct genomes
	hits <- blast_reciprocal_hits(blast, blast, min_cover = min_cover) |>
		filter(str_remove(query_id, ":[^:]*$") != str_remove(subject_id, ":[^:]*$"))
  
	# Clusters of genes linked by their reciprocal hits
	g <- tbl_graph(
		nodes = genes,
		edges = hits,
		directed = FALSE,
		node_key = "gene_uid"
	) |>
		mutate(cluster_id = group_components())
		
	# Keep a single paralog per genome, the one with best centrality
	g <- g |> 
		mutate(score = centrality_closeness(weights = pmax(1e-6,100 - pct_ident))) |>
		group_by(cluster_id, genome) |>
		mutate(rk = row_number(score)) |> ungroup() |> filter(rk<=1L) |> select(!rk) |> # equivalent to filter(row_number(score)<=1L)
		ungroup()
	
	# Determine template genes: one per core cluster and with the best closeness score
	g <- g |> 
		group_by(cluster_id) |>
		mutate(is_core_template = (n_distinct(gene_id)==n_distinct(genome)) & (n_distinct(genome)>=n_core) & (row_number(score)<=1L)) |>
		ungroup()

	#ggraph(g) + geom_edge_link(color="grey") + geom_node_point(aes(color=is_core_template)) + scale_color_manual(values = c("TRUE"="black","FALSE"="grey")) + coord_equal()
  #g %>% filter(is_core_template)
  
  return(g)
}

