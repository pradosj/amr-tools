# Functions of FASTA_ROTATE (sourced after lib_io.R, whose read_blast6() and read_fasta() it uses)

suppressPackageStartupMessages({
	library(GenomicRanges)
	library(Biostrings)
	library(tidyverse)
})


seq_rotate <- function(fa,rotation_pos) {
  FA <- xscat(subseq(fa, start = start(rotation_pos)), subseq(fa, end = end(rotation_pos)))
  minus <- as.logical(strand(rotation_pos) == "-")
  FA[minus] <- reverseComplement(FA[minus])
  names(FA) <- names(fa)
  FA
}

#' Rotate the sequences to start at their best hit of the reference
#'
#' Each sequence starts at its best hit (highest bit score) of the reference, on the strand of the hit;
#' the sequences without hit are not rotated.
#'
#' @param qry_fa FASTA file of the sequences to rotate.
#' @param blast Hits of the sequences to rotate (query) on the reference (subject), as read by read_blast6().
fasta_rotate_blast <- function(qry_fa, blast) {
	fa <- read_fasta(qry_fa)
	hits <- read_blast6(blast) |>
		group_by(query_id) |>
		slice_max(bit_score, n = 1, with_ties = FALSE) |>    # ties: the first hit, as sorted by BLAST
		ungroup()
	rotation_pos <- GRanges(
		hits$query_id,
		IRanges(hits$query_start, hits$query_end),    # BLAST reports the query on its + strand
		strand = ifelse(hits$subject_start > hits$subject_end, "-", "+")
	)
	rotation_pos <- c(rotation_pos, GRanges(names(fa), IRanges(1, width = 1), strand = "+"))
	rotation_pos <- rotation_pos[match(names(fa), as.character(seqnames(rotation_pos)))]
	rotation_pos <- resize(rotation_pos, width = 0)

	seq_rotate(fa,rotation_pos)  
}



fasta_rotate_bam <- function(qry_fa,bam) {
    aln <- readGAlignments(bam,use.names=TRUE)
    
    # For each query, determine the ref it cover most
    qry_map <- split(aln,names(aln)) |> 
      lapply(\(e) sum(coverage(e)>0)) |>
      enframe(name="query_id") |>
      unnest_longer(value,indices_to="ref_id",values_to="coverage") |>
      group_by(query_id) |>
      slice_max(coverage,n=1,with_ties=FALSE)
    
    # Filter BAM to keep best contig pairs, and longest alignments
    aln <- subset(aln,paste(names(aln),seqnames(aln)) %in% paste(qry_map$query_id,qry_map$ref_id))
    aln <- aln[order(width(aln),decreasing=TRUE)]
    aln <- aln[!duplicated(names(aln))]
    
    # Use selected alignment to rotate queries
    rotation_pos <- range(bam_cigar_to_bedpe(aln)$query) |> resize(width=0L)
    
    fa <- read_fasta(qry_fa)
    seq_rotate(fa,rotation_pos)
}
