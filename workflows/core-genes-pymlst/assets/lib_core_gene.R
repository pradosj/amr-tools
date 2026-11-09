
source("assets/lib_io.R")


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
core_genes <- function(hits_tbl, ffn_tbl, n_core) {
  library(tidygraph)
  
  # Read ffn_tbl to determine all gene ids, not only the ones in hits
	genes <- ffn_tbl %>% 
		mutate(gene_id = map(file,~{
			seq <- read_fasta(.x)
			names(seq)
		})) %>% 
		select(!file) %>% 
		unnest(-1) %>% 
		mutate(gene_uid = str_c(genome,":",gene_id)) %>% 
		relocate(gene_uid)
  
  # Read hits
  hits <- hits_tbl |>
    mutate(blast = map(file, read_blast6)) |>
    unnest(blast) |>
  	left_join(select(genes,query_uid=gene_uid,query_id=gene_id,query_genome=genome),by=c("query_genome","query_id")) %>% 
  	left_join(select(genes,subject_uid=gene_uid,subject_id=gene_id,subject_genome=genome),by=c("subject_genome","subject_id")) %>%
    relocate(query_uid,subject_uid)
  
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

