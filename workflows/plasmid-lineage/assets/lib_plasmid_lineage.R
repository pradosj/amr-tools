# Functions of the plasmid-lineage pipeline, used by its R processes and by its report (plasmid_lineage.qmd)

suppressPackageStartupMessages({
	library(GenomicRanges)
	library(Biostrings)
	library(tidyverse)
})


# ---------------------------------------
#  Input / output
# ---------------------------------------

# BLAST tabular output, as produced with -outfmt "6 std qlen slen"
read_blast6 <- function(f) {
	readr::read_tsv(
		f,
		col_names = c(
			"query_id","subject_id","pct_ident","align_len","num_mismatch",
			"num_gapopen","query_start","query_end","subject_start","subject_end",
			"evalue","bit_score","query_len","subject_len"
		),
		col_types = "ccdiiiiiiiddii"
	)
}

# Read a FASTA file, naming the sequences by their id (first word of the header) as BLAST does
read_fasta <- function(f) {
	fa <- readDNAStringSet(f)
	names(fa) <- sub("\\s.*", "", names(fa))
	fa
}

# File name of a sequence id
seq_file_name <- function(seq_id, ext = ".fasta") {
	paste0(make.names(seq_id), ext)
}

# Minimum number of plasmids to call a gene or a region core: the given count, or 80% of the plasmids
core_count <- function(n_core, n_plasmids) {
	if (is.null(n_core) || is.na(n_core)) ceiling(0.8 * n_plasmids) else as.integer(n_core)
}


# ---------------------------------------
#  Alignments
# ---------------------------------------

subject_granges <- function(x, ...) {
	GRanges(
		x$subject_id,
		IRanges(pmin(x$subject_start, x$subject_end), pmax(x$subject_start, x$subject_end)),
		strand = ifelse(x$subject_start > x$subject_end, "-", "+"),
		...
	)
}

query_granges <- function(x, ...) {
	GRanges(
		x$query_id,
		IRanges(x$query_start, x$query_end),
		strand = "+",
		...
	)
}

# Regions of the queries, with the subjects aligning on each of them
find_shared_regions <- function(aln) {
	d <- query_granges(aln) |>
		unstrand() |>
		disjoin()
	h <- findOverlaps(d, unstrand(query_granges(aln)), type = "within")
	d$subjects <- extractList(aln$subject_id, as(h, "List")) |> unique()
	d
}

# Remove the alignments of a query included in another alignment of the same query on the same subject region
filter_inner_alignments <- function(aln) {
	S <- subject_granges(aln) |> unstrand()
	w <- findOverlaps(S, S, type = "within")
	w <- w[queryHits(w) != subjectHits(w)]
	w <- w[aln$query_id[queryHits(w)] == aln$query_id[subjectHits(w)]]
	filter(aln, !seq(nrow(aln)) %in% queryHits(w))
}

# Polygons joining the aligned regions of the query and of the subject
aln_polygons <- function(aln) {
	poly <- aln |>
		mutate(poly_id = seq(n()))
	bind_rows(
		poly |>
			select(poly_id, seq = query_id, query_start, query_end, pct_ident, num_mismatch, num_gapopen) |>
			pivot_longer(c(query_start, query_end)),
		poly |>
			select(poly_id, seq = subject_id, subject_start, subject_end, pct_ident, num_mismatch, num_gapopen) |>
			pivot_longer(c(subject_start, subject_end))
	) |>
		mutate(name = fct_relevel(name, "query_start", "query_end", "subject_end", "subject_start")) |>
		arrange(poly_id, name)
}

# Mean identity of the alignments of each pair of sequences
blast_ani <- function(aln, min_len = 500) {
	aln |>
		filter(align_len >= min_len) |>
		group_by(query_id, subject_id) |>
		summarize(ani = mean(pct_ident), .groups = "drop")
}

# Clustering of the sequences on their ANI
ani_hclust <- function(ani) {
	M <- ani |>
		pivot_wider(id_cols = query_id, names_from = subject_id, values_from = ani, values_fill = 0) |>
		column_to_rownames("query_id") |>
		as.matrix()
	M <- M[, rownames(M), drop = FALSE]
	as.dist(100 - M) |>
		hclust(method = "single")
}


# ---------------------------------------
#  Plots
# ---------------------------------------

# Dot plot of the alignments of each pair of sequences
gg_pw_alignment <- function(aln) {
	aln |>
		ggplot(aes(x = query_start, xend = query_end, y = subject_start, yend = subject_end, color = pct_ident)) +
			facet_grid(subject_id ~ query_id) +
			geom_segment() +
			coord_equal() +
			theme(panel.spacing = unit(0, "mm")) +
			theme(strip.text.y.right = element_text(angle = 0))
}

# Polygons of the alignments between consecutive sequences, in the order of `seq_levels`
alignment_polygons <- function(aln, seq_levels, min_len = 1000) {
	aln |>
		mutate(across(c(query_id, subject_id), ~ factor(., seq_levels))) |>
		filter(abs(as.integer(query_id) - as.integer(subject_id)) == 1L) |>
		filter(align_len >= min_len) |>
		filter_inner_alignments() |>
		aln_polygons() |>
		mutate(seq = factor(seq, seq_levels)) |>
		group_by(poly_id) |>
		filter(max(as.integer(seq)) - min(as.integer(seq)) == 1L) |>
		ungroup()
}

# Plot of the polygons of alignment_polygons(), with the sequences as horizontal lines
gg_alignment_polygons <- function(polygons, seqlen) {
	centers <- polygons |>
		group_by(poly_id) |>
		summarise(x = mean(value), y = mean(as.numeric(seq)), label = str_c(first(num_mismatch), "|", first(num_gapopen)))
	ggplot(polygons) +
		geom_polygon(aes(x = value, y = seq, group = poly_id, fill = pct_ident), alpha = 0.66, color = "black") +
		geom_text(aes(x = x, y = y, label = label), data = centers, size = 3) +
		geom_segment(aes(x = 0, xend = seqlen, y = seq, yend = seq), data = seqlen, linewidth = 2) +
		scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
		scale_fill_continuous(limits = c(98, 100), oob = scales::squish, low = "red") +
		theme_bw() +
		ylab("") + xlab("Genomic position")
}


# ---------------------------------------
#  Steps of the pipeline
# ---------------------------------------

# Longest region of a sequence shared by at least min_frac of the sequences, written to `out`
# (empty when there is none): the rotated sequences start at this region
find_align_seq <- function(fasta, blast, out, min_len = 1000, min_frac = 0.8) {
	fa <- read_fasta(fasta)
	aln <- read_blast6(blast) |>
		filter(align_len >= min_len)
	align_seq <- DNAStringSet()
	if (nrow(aln) > 0) {
		d <- find_shared_regions(aln)
		d <- d[lengths(d$subjects) / n_distinct(aln$subject_id) >= min_frac]
		if (length(d) > 0) align_seq <- fa[d[which.max(width(d))]]
	}
	writeXStringSet(align_seq, out)
}

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

# Core genome: regions of each plasmid shared by at least n_core plasmids (itself included). Writes:
#   - core_genome.tsv: the regions of each plasmid, with the number of plasmids sharing it and whether it is core
#   - fastani_genomes/<plasmid>.fasta: the regions of each plasmid, for FastANI
#   - fastani_genomes.tsv: the plasmid of each FastANI genome
core_genome <- function(fasta, blast, n_core = NA, min_len = 500) {
	FA <- read_fasta(fasta)
	n_core <- core_count(n_core, length(FA))
	hits <- read_blast6(blast) |>
		filter(align_len >= min_len) |>
		filter(query_id != subject_id)
	hits <- c(
		GRanges(hits$query_id, IRanges(hits$query_start, hits$query_end), name = hits$subject_id, seqinfo = Seqinfo(names(FA), width(FA))),
		GRanges(hits$subject_id, IRanges(pmin(hits$subject_start, hits$subject_end), pmax(hits$subject_start, hits$subject_end)), name = hits$query_id, seqinfo = Seqinfo(names(FA), width(FA)))
	)
	hits <- splitAsList(hits, hits$name) |>
		GenomicRanges::reduce() |>
		coverage() |>
		as("GRanges")
	hits$n_plasmids <- hits$score + 1L    # the plasmids aligning, plus itself
	core_mask <- hits |>
		subset(n_plasmids >= n_core) |>
		GenomicRanges::reduce() |>
		subset(width >= min_len) |>
		gaps(ignore.strand = TRUE)
	hits$is_core <- countOverlaps(hits, core_mask) == 0
	hits |>
		as_tibble() |>
		select(seqnames, start, end, n_plasmids, is_core) |>
		write_tsv("core_genome.tsv")

	genomes <- tibble(plasmid_id = seqlevels(hits), file = seq_file_name(plasmid_id))
	dir.create("fastani_genomes")
	for (i in seq_len(nrow(genomes))) {
		h <- hits[seqnames(hits) == genomes$plasmid_id[i]]
		dna <- FA[h]
		names(dna) <- as.character(h)
		writeXStringSet(dna, file.path("fastani_genomes", genomes$file[i]))
	}
	write_tsv(genomes, "fastani_genomes.tsv")
}
