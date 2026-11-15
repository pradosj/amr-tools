
library(GenomicAlignments)
library(tidyverse)


bam_cigar_to_bedpe <- function(aln) {
	stopifnot("aln must be named"=!is.null(names(aln)))
	
	ref_rg <- cigarRangesAlongReferenceSpace(
		cigar(aln), pos = start(aln),
		ops = c("M", "=", "X"), reduce.ranges = FALSE,f=as.factor(seqnames(aln))
	) |>
		GRanges(seqinfo = seqinfo(aln))
	
	
	qry_len <- cigarWidthAlongQuerySpace(cigar(aln), before.hard.clipping = TRUE)
	qry_rg <- cigarRangesAlongQuerySpace(
		cigar(aln),
		ops = c("M", "=", "X"), reduce.ranges = FALSE,
		before.hard.clipping = TRUE
	) |>
		setNames(names(aln))
	aln_idx <- rep(seq_along(aln),lengths(qry_rg))
	
	qry_rg <- GRanges(
		qry_rg,
		strand = rep(strand(aln),lengths(qry_rg)),
		seqinfo = data.frame(seqnames=names(aln),seqlengths=qry_len,stringsAsFactors = FALSE) |> distinct() |> as("Seqinfo")
	)
	ranges(qry_rg)[strand(qry_rg)=="-"] <- reflect(
		ranges(qry_rg)[strand(qry_rg)=="-"],
		IRanges(1,seqlengths(qry_rg)[as.integer(seqnames(qry_rg[strand(qry_rg)=="-"]))])
	)
	
	DataFrame(aln_idx,ref = ref_rg,query = qry_rg)
}



# Generate one polygon per alignment
bedpe_polygons <- function(bedpe) {
	bind_rows(
		"ref" = as_tibble(bedpe$ref) |> bind_cols(aln_idx = bedpe$aln_idx, bedpe_idx=seq_along(bedpe$aln_idx)),
		"qry" = as_tibble(bedpe$query) |> bind_cols(aln_idx = bedpe$aln_idx, bedpe_idx=seq_along(bedpe$aln_idx)),
		.id = "source"
	) |>
		select(-width) |>
		mutate(source = factor(source,c("ref","qry"))) |>
		pivot_longer(cols = c("start","end"),names_to = "order",names_transform = list(order=~as.integer(factor(.,c("start","end"))))) |>
		mutate(order = case_when(
			source=="ref" ~ order,
			source=="qry" & strand!="-" ~ 5 - order,
			source=="qry" & strand=="-" ~ 2 + order,
			TRUE ~ NA_integer_
		)) |>
		select(-strand) |>
		arrange(aln_idx,bedpe_idx,order)
}





