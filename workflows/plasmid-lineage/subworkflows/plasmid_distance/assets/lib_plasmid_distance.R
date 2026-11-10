

filter_inner_alignments <- function(ALN) {
	S <- subject_granges(ALN) |> unstrand()
	#Q <- query_granges(ALN) |> unstrand()
	# Find alignments coming from the same query sequence and targeting the same subject region
	w <- findOverlaps(S,S,type="within")
	w <- w[queryHits(w)!=subjectHits(w)]
	w <- w[ALN$query_id[queryHits(w)] == ALN$query_id[subjectHits(w)]]
	
	filter(ALN,!seq(nrow(ALN)) %in% queryHits(w))
}


blast_distance_graph <- function(blast, min_len = 1000) {
		ALN <- read_blast6(blast)
		
		g <- ALN |>
			filter(align_len>=min_len) |>
			filter_inner_alignments() |>
			group_by(query_id,subject_id) |>
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
				qid = pmin(as.character(query_id),as.character(subject_id)),
				sid = pmax(as.character(query_id),as.character(subject_id))
			) |>
			slice_max(weight,n = 1L,with_ties = FALSE)
		
		tbl_graph(
			nodes = tibble(name=unique(c(g$qid,g$sid))),
			edges = g,
			directed = FALSE
		)
}

