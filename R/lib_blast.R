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

#' Ranges of the subjects of alignments, on the strand of the alignment
#'
#' @param x Alignments, as read by read_blast6().
#' @param ... Metadata columns of the GRanges.
#' @return A GRanges: one range per alignment, on its subject.
subject_granges <- function(x, ...) {
	GenomicRanges::GRanges(
		x$subject_id,
		IRanges::IRanges(pmin(x$subject_start, x$subject_end), pmax(x$subject_start, x$subject_end)),
		strand = ifelse(x$subject_start > x$subject_end, "-", "+"),
		...
	)
}

#' Ranges of the queries of alignments (strand +)
#'
#' @param x Alignments, as read by read_blast6().
#' @param ... Metadata columns of the GRanges.
#' @return A GRanges: one range per alignment, on its query.
query_granges <- function(x, ...) {
	GenomicRanges::GRanges(
		x$query_id,
		IRanges::IRanges(x$query_start, x$query_end),
		strand = "+",
		...
	)
}

#' Remove the alignments included in another alignment of the same query on the same subject region
#'
#' @param aln Alignments, as read by read_blast6().
#' @return The alignments that are not included in another one.
blast_filter_inner_alignments <- function(aln) {
	S <- BiocGenerics::unstrand(subject_granges(aln))
	w <- IRanges::findOverlaps(S, S, type = "within")
	w <- w[S4Vectors::queryHits(w) != S4Vectors::subjectHits(w)]
	w <- w[aln$query_id[S4Vectors::queryHits(w)] == aln$query_id[S4Vectors::subjectHits(w)]]
	filter(aln, !seq(nrow(aln)) %in% S4Vectors::queryHits(w))
}

#' Distance graph of the sequences, from their BLASTN against themselves
#'
#' For each pair of sequences, the alignments of at least min_len bp that are not included in another one (see
#' filter_inner_alignments()) give the errors: the mismatches, plus one per extra fragment (n_frag - 1). The weight of
#' an edge is the number of errors per aligned bp, the worst of the two directions.
#'
#' @param aln Blast alignments, as read by read_blast6().
#' @return An undirected tbl_graph: the sequences as nodes (name), and an edge per pair of sequences aligning (each
#'   sequence with itself included), with num_mis, num_gapopen, n_frag, total_frag_len, err and weight.
blast_distance_graph <- function(aln) {
	g <- aln |>
		group_by(query_id, subject_id) |>
		summarise(
			num_mis = sum(num_mismatch),
			num_gapopen = sum(num_gapopen),
			n_frag = n(),
			total_frag_len = sum(align_len)
		) |>
		mutate(
			err = (num_mis + (n_frag - 1)),
			weight = err / total_frag_len
		) |>
		group_by(
			qid = pmin(as.character(query_id), as.character(subject_id)),
			sid = pmax(as.character(query_id), as.character(subject_id))
		) |>
		slice_max(weight, n = 1L, with_ties = FALSE)

	tbl_graph(
		nodes = tibble(name = unique(c(g$qid, g$sid))),
		edges = g,
		directed = FALSE
	)
}



