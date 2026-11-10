# Functions of the plasmid-lineage pipeline, used by its R processes and by its report (plasmid_lineage.qmd)

suppressPackageStartupMessages({
	library(GenomicRanges)
	library(Biostrings)
	library(tidyverse)
})

# ---------------------------------------
#  Input / output
# ---------------------------------------

# ---------------------------------------
#  Alignments
# ---------------------------------------

#' Regions of the queries, with the subjects aligning on each of them
#'
#' @param aln Alignments, as read by read_blast6().
#' @return A GRanges of the disjoint regions of the queries, with the ids of the subjects aligning on each of
#'   them in `subjects`.
find_shared_regions <- function(aln) {
	d <- query_granges(aln) |>
		unstrand() |>
		disjoin()
	h <- findOverlaps(d, unstrand(query_granges(aln)), type = "within")
	d$subjects <- extractList(aln$subject_id, as(h, "List")) |> unique()
	d
}

#' Polygons joining the aligned regions of the query and of the subject
#'
#' @param aln Alignments, as read by read_blast6().
#' @return A tibble of the polygon corners: poly_id, seq, name (corner), value (position).
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

#' Mean identity of the alignments of each pair of sequences
#'
#' @param aln Alignments, as read by read_blast6().
#' @param min_len Minimum length of the alignments.
#' @return A tibble with query_id, subject_id and ani.
blast_ani <- function(aln, min_len = 500) {
	aln |>
		filter(align_len >= min_len) |>
		group_by(query_id, subject_id) |>
		summarize(ani = mean(pct_ident), .groups = "drop")
}

#' Single linkage clustering of the sequences on their ANI
#'
#' @param ani ANI of the pairs of sequences, as computed by blast_ani().
#' @return An hclust object, on the distance 100 - ANI.
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

#' Dot plot of the alignments of each pair of sequences
#'
#' @param aln Alignments, as read by read_blast6().
#' @return A ggplot.
gg_pw_alignment <- function(aln) {
	aln |>
		ggplot(aes(x = query_start, xend = query_end, y = subject_start, yend = subject_end, color = pct_ident)) +
			facet_grid(subject_id ~ query_id) +
			geom_segment() +
			coord_equal() +
			theme(panel.spacing = unit(0, "mm")) +
			theme(strip.text.y.right = element_text(angle = 0))
}

#' Polygons of the alignments between consecutive sequences
#'
#' @param aln Alignments, as read by read_blast6().
#' @param seq_levels The sequences, in their plotting order.
#' @param min_len Minimum length of the alignments.
#' @return The polygons, as computed by aln_polygons(), of the alignments between consecutive sequences of seq_levels.
alignment_polygons <- function(aln, seq_levels, min_len = 1000) {
	aln |>
		mutate(across(c(query_id, subject_id), ~ factor(., seq_levels))) |>
		filter(abs(as.integer(query_id) - as.integer(subject_id)) == 1L) |>
		filter(align_len >= min_len) |>
		blast_filter_inner_alignments() |>
		aln_polygons() |>
		mutate(seq = factor(seq, seq_levels)) |>
		group_by(poly_id) |>
		filter(max(as.integer(seq)) - min(as.integer(seq)) == 1L) |>
		ungroup()
}

#' Plot of the alignment polygons, with the sequences as horizontal lines
#'
#' @param polygons Polygons, as computed by alignment_polygons().
#' @param seqlen A tibble with the sequences (seq) and their length (seqlen).
#' @return A ggplot.
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

#' Longest region of a sequence shared by most of the sequences
#'
#' The pipeline rotates the sequences to start at this region.
#'
#' @param fasta FASTA file of the sequences.
#' @param blast Self alignment of the sequences, as read by read_blast6().
#' @param out FASTA file written with the region (empty when there is none).
#' @param min_len Minimum length of the alignments.
#' @param min_frac Minimum fraction of the sequences sharing the region.
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

