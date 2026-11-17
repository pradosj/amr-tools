# Functions of FIND_ROTATION_ORIGINS (sourced after lib_io.R and lib_blast.R)

suppressPackageStartupMessages({
	library(GenomicRanges)
	library(Biostrings)
	library(tidyverse)
})


#' Regions of the queries, with the subjects aligning on each of them
#'
#' @param aln Alignments, as read by read_blast6().
#' @return A GRanges of the disjoint regions of the queries, with the ids of the subjects aligning on each of
#'   them in `subjects`.
find_shared_regions <- function(aln) {
	d <- blast_query_granges(aln) |>
		unstrand() |>
		disjoin()
	h <- findOverlaps(d, unstrand(blast_query_granges(aln)), type = "within")
	d$subjects <- extractList(aln$subject_id, as(h, "List")) |> unique()
	d
}

#' Origins of the rotation of a set of genomes
#'
#' The regions shared by most of the genomes, chosen from the longest one: a region is kept when none of the
#' sequences aligning on it already aligns on a kept region, so that each replicon (chromosome, plasmid) of the
#' genomes gets at most one origin.
#'
#' @param fasta FASTA file of the sequences of all the genomes, named `<genome id>:<sequence id>` (FASTA_CAT).
#' @param blast Self alignment of the sequences, as read by read_blast6().
#' @param out FASTA file written with the origins (empty when there is none).
#' @param min_len Minimum length of the alignments.
#' @param min_frac Minimum fraction of the genomes sharing a region.
#' @return The origins (GRanges on the sequences of fasta), invisibly.
find_rotation_origins <- function(fasta, blast, out, min_len = 1000, min_frac = 0.8) {
	genome_of <- function(seq_id) str_replace(seq_id, ":[^:]*$", "")
	fa <- read_fasta(fasta)
	names(fa) <- str_replace(names(fa), " .*", "")
	n_genomes <- n_distinct(genome_of(names(fa)))
	aln <- read_blast6(blast) |>
		filter(align_len >= min_len)
	origins <- GRanges()
	if (nrow(aln) > 0) {
		d <- find_shared_regions(aln)
		genomes <- relist(genome_of(unlist(d$subjects)), d$subjects) |> unique()
		d <- d[lengths(genomes) / n_genomes >= min_frac]
		d <- d[order(width(d), decreasing = TRUE)]
		covered <- character(0)
		keep <- logical(length(d))
		for (i in seq_along(d)) {
			s <- d$subjects[[i]]
			if (!any(s %in% covered)) {
				keep[i] <- TRUE
				covered <- c(covered, s)
			}
		}
		origins <- d[keep]
	}
	origins_fa <- if (length(origins) > 0) fa[origins] else DNAStringSet()
	writeXStringSet(origins_fa, out)
	invisible(origins)
}
