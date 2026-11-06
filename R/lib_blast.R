# Functions on BLAST results, shared by the R steps of the pipelines (reached through file links; sourced after
# lib_io.R, whose read_blast6() reads the BLAST results they take)

#' Reciprocal BLAST hits of two sets of sequences A and B
#'
#' The pairs of sequences of A and B aligned in both directions (hits covering more than min_cover of the query
#' sequence), with the identity and score of the worst direction (lowest bit score) of each pair.
#'
#' @param blast1 BLAST (BLASTP or BLASTN) of A (query) against B (subject), in the format read by read_blast6().
#' @param blast2 BLAST of B (query) against A (subject), in the format read by read_blast6(). With A = B (blast1 =
#'   blast2), each pair of distinct sequences is given in both orders, and each sequence is paired with itself.
#' @param min_cover Minimum coverage of the query sequence by its alignment.
#' @return The hits of blast1 aligned in both directions, one per pair (query in A, subject in B), in the format of
#'   read_blast6(): pct_ident, align_len, num_mismatch, num_gapopen, evalue and bit_score of the worst direction; the
#'   positions and lengths of blast1.
blast_reciprocal_hits <- function(blast1, blast2, min_cover = 0.75) {
	# Best hit of each pair of query and subject
	best_hits <- function(blast) {
		read_blast6(blast) |>
			filter(align_len / query_len > min_cover) |>
			group_by(query_id, subject_id) |>
			slice_max(bit_score, n = 1, with_ties = FALSE) |>
			ungroup()
	}
	hits1 <- best_hits(blast1)
	hits2 <- best_hits(blast2)

	# Pairs aligned in both directions, with the worst direction
	inner_join(hits1, hits2, by = c(query_id = "subject_id", subject_id = "query_id"), suffix = c("", "_2")) |>
		mutate(
			worst_2 = bit_score_2 < bit_score,    # the second direction (blast2) is the worst
			# the scores of the worst direction
			across(
				c(pct_ident, align_len, num_mismatch, num_gapopen, evalue, bit_score),
				\(x) if_else(worst_2, get(str_c(cur_column(), "_2")), x)
			)
		) |>
		select(!ends_with("_2"))
}
