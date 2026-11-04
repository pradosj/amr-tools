# Rotation of genome sequences on their hits of a reference dataset (ROTATE_GENOMES)

suppressPackageStartupMessages({
	library(GenomicRanges)
	library(Biostrings)
	library(tidyverse)
})

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

# Rotate the sequences to start at their best hit (highest bit score) of the reference, on the strand of
# the hit (sequences without hit are not rotated). `blast`: hits of the reference sequences (query) on the
# sequences to rotate (subject).
rotate_fasta <- function(fasta, blast, out) {
	fa <- read_fasta(fasta)
	hits <- read_blast6(blast) |>
		group_by(subject_id) |>
		slice_max(bit_score, n = 1, with_ties = FALSE) |>    # ties: the first hit, as sorted by BLAST
		ungroup()
	rotation_pos <- GRanges(
		hits$subject_id,
		IRanges(pmin(hits$subject_start, hits$subject_end), pmax(hits$subject_start, hits$subject_end)),
		strand = ifelse(hits$subject_start > hits$subject_end, "-", "+")
	)
	rotation_pos <- c(rotation_pos, GRanges(names(fa), IRanges(1, width = 1), strand = "+"))
	rotation_pos <- rotation_pos[match(names(fa), as.character(seqnames(rotation_pos)))]
	rotation_pos <- resize(rotation_pos, width = 0)

	FA <- xscat(subseq(fa, start = start(rotation_pos)), subseq(fa, end = end(rotation_pos)))
	minus <- as.logical(strand(rotation_pos) == "-")
	FA[minus] <- reverseComplement(FA[minus])
	names(FA) <- names(fa)
	writeXStringSet(FA, out)
}
