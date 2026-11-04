# Functions of ROTATE_GENOMES (sourced after lib_plasmid_lineage.R, whose functions it uses)

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
