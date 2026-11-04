# Functions of the CORE_GENES_PYMLST subworkflow (sourced after lib_plasmid_lineage.R, whose functions it uses)

# Core genes: clusters of reciprocal best BLASTP hits of the prokka proteins, with one gene in at least
# n_core plasmids. Writes:
#   - core_genes.tsv: the genes of the clusters, with their cluster and whether the cluster is core
#   - core_clusters.tsv: the size and number of plasmids of each cluster
#   - pymlst_genomes/<plasmid>.fasta: the core genes of each plasmid, for pyMLST
#   - pymlst_template.fasta: the core genes of the plasmid with the most core genes, template of the pyMLST database
#   - pymlst_strains.tsv: the plasmid of each pyMLST strain
core_genes <- function(prokka_dir, blastp, n_core = NA, min_cover = 0.75) {
	library(tidygraph)
	gff <- rtracklayer::import.gff3(file.path(prokka_dir, "prokka.gff"))
	gff <- gff[gff$type == "CDS"]
	n_core <- core_count(n_core, length(seqlevelsInUse(gff)))

	# Reciprocal hits, keeping the worst direction of each pair of genes
	z <- read_blast6(blastp) |>
		filter(align_len / query_len > min_cover) |>
		group_by(query_id, subject_id) |>
		slice_max(bit_score, n = 1, with_ties = FALSE) |>
		ungroup() |>
		filter(str_c(query_id, subject_id, sep = ":") %in% str_c(subject_id, query_id, sep = ":")) |>
		filter(query_id != subject_id) |>
		mutate(undirected_edge_id = str_c(pmin(query_id, subject_id), pmax(query_id, subject_id), sep = ":")) |>
		group_by(undirected_edge_id) |>
		slice_min(bit_score, n = 1, with_ties = FALSE) |>
		ungroup()

	# Clusters of genes, keeping one paralog per plasmid in each cluster
	g <- tbl_graph(
		nodes = tibble(ID = gff$ID, contig_id = as.character(seqnames(gff)), width = width(gff)),
		edges = select(z, from = query_id, to = subject_id, bit_score, pct_ident),
		directed = FALSE,
		node_key = "ID"
	) |>
		mutate(cluster_id = group_components()) |>
		mutate(score = centrality_closeness(weights = 100.1 - pct_ident)) |>
		group_by(cluster_id, contig_id) |>
		mutate(rk = row_number(score)) |>
		ungroup() |>
		filter(rk <= 1)

	genes <- as_tibble(g)
	clusters <- genes |>
		group_by(cluster_id) |>
		summarize(cluster_size = n_distinct(ID), distinct_plasmid = n_distinct(contig_id), mean_width = mean(width)) |>
		filter(cluster_size == distinct_plasmid) |>
		mutate(is_core = distinct_plasmid >= n_core)
	genes <- genes |>
		select(ID, contig_id, width, cluster_id) |>
		left_join(select(clusters, cluster_id, is_core), by = "cluster_id") |>
		mutate(is_core = replace_na(is_core, FALSE))
	write_tsv(genes, "core_genes.tsv")
	write_tsv(clusters, "core_clusters.tsv")

	# Sequences of the core genes of each plasmid
	ffn <- read_fasta(file.path(prokka_dir, "prokka.ffn"))
	core <- filter(genes, is_core)
	dna <- split(ffn[core$ID], core$contig_id)
	strains <- tibble(contig_id = names(dna), file = seq_file_name(contig_id)) |>
		mutate(strain = fs::path_ext_remove(file))
	dir.create("pymlst_genomes")
	for (i in seq_len(nrow(strains))) writeXStringSet(dna[[strains$contig_id[i]]], file.path("pymlst_genomes", strains$file[i]))
	template <- names(dna)[which.max(lengths(dna))]
	writeXStringSet(if (length(template)) dna[[template]] else DNAStringSet(), "pymlst_template.fasta")
	write_tsv(select(strains, strain, contig_id), "pymlst_strains.tsv")
}
