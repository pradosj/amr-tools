# Functions on BLAST results, shared by the R steps of the pipelines (reached through file links; sourced after
# lib_io.R, whose read_blast6() reads the BLAST results they take)

#' Reciprocal BLASTP hits of two sets of sequences A and B
#'
#' The pairs of sequences of A and B aligned in both directions (hits covering more than min_cover of the query
#' protein), with the identity and score of the worst direction (lowest bit score) of each pair.
#'
#' @param blast1 BLASTP of A (query) against B (subject), in the format read by read_blast6().
#' @param blast2 BLASTP of B (query) against A (subject), in the format read by read_blast6(). With A = B (blast1 =
#'   blast2), each pair of distinct sequences is given in both orders, and each sequence is paired with itself.
#' @param min_cover Minimum coverage of the query protein by its alignment.
#' @return A tibble with one row per pair: query_id (in A), subject_id (in B), and the pct_ident, align_len, evalue and
#'   bit_score of the worst direction.
blast_reciprocal_hits <- function(blast1, blast2, min_cover = 0.75) {
	# Best hit of each pair of query and subject
	best_hits <- function(blast) {
		read_blast6(blast) |>
			filter(align_len / query_len > min_cover) |>
			group_by(query_id, subject_id) |>
			slice_max(bit_score, n = 1, with_ties = FALSE) |>
			ungroup() |>
			select(query_id, subject_id, pct_ident, align_len, evalue, bit_score)
	}
	hits1 <- best_hits(blast1)
	hits2 <- best_hits(blast2)

	# Pairs aligned in both directions, with the worst direction
	inner_join(hits1, hits2, by = c(query_id = "subject_id", subject_id = "query_id"), suffix = c("", "_2")) |>
		mutate(
			worst_2   = bit_score_2 < bit_score,
			pct_ident = if_else(worst_2, pct_ident_2, pct_ident),
			align_len = if_else(worst_2, align_len_2, align_len),
			evalue    = if_else(worst_2, evalue_2, evalue),
			bit_score = if_else(worst_2, bit_score_2, bit_score)
		) |>
		select(!ends_with("_2"))
}
