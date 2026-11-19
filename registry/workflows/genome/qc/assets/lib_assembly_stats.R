# QC statistics of an assembly (COMBINE_QC_STATS, written as qc.json), its dashboard (genome_qc_report.qmd)
# and the table of all the assemblies (GENOME_MULTIQC). Sourced after lib_io.R.

library(tidyverse)
library(Biostrings)

# Extract statistics of each contig from the sequence
fa_extract_contigs_stats <- function(fa_file) {
	fa <- readDNAStringSet(fa_file)
	tibble(
		contig_idx = seq_along(fa),
		contig_id = str_remove(names(fa)," .*"),
		contig_name = names(fa),
		length = lengths(fa),
		GC_count = as.vector(Biostrings::letterFrequency(fa,"GC")),
		N_count = as.vector(Biostrings::letterFrequency(fa,"N"))
	)
}

#' Summary of the alignment of a set of reads on an assembly
samtools_stats_summary <- function(stats_file) {
	s <- read_samtools_stats(stats_file)
	enframe(s) %>% 
	  pivot_wider() %>% 
	  mutate(
	    num_sequenced = as.numeric(pluck(SN, "raw.total.sequences", .default = NA_real_)),
	    num_mapped = as.numeric(pluck(SN, "reads.mapped", .default = NA_real_)),
	    num_mapped_bp = as.numeric(pluck(SN, "bases.mapped..cigar.", .default = NA_real_))
	  ) %>% 
	  mutate(
	    pct_mapped = 100 * num_mapped / num_sequenced,
	    num_mapped_bp = num_mapped_bp,
	    depth_hist = map(COV_HIST,~as_tibble(.) |> select(depth=range_max,count))
	  )
}

#' Merge all statistics of an assembly found in dir
merge_genome_qc_stats <- function(dir) {
  
  # Get basic contigs statistics from the sequence
	contigs <- fa_extract_contigs_stats(fs::path(dir,'assembly.fasta'))
	
	# Add read coverage statistics to the contigs
	contigs <- local({
	  long_cov <- read_samtools_coverage(fs::path(dir,'long_reads.coverage.tsv'))
	  short_cov <- read_samtools_coverage(fs::path(dir,'short_reads.coverage.tsv'))
	  contigs |>
	    left_join(select(long_cov,contig_id,long_meandepth = meandepth),by = "contig_id") |>
	    left_join(select(short_cov,contig_id,short_meandepth = meandepth),by = "contig_id")
	})

	# Compute summary assembly statistics from contigs
	assembly <- contigs %>% 
	  summarize(
	    num_contig = n(),
	    total_len = sum(length),
	    min_len = min(length),
	    median_len = median(length),
	    max_len = max(length),
	    N50 = N50(length),
	    GC_pct = 100 * sum(GC_count) / total_len,
	    long_meandepth = weighted.mean(long_meandepth,length),
	    short_meandepth = weighted.mean(short_meandepth,length)
    )

	# Add checkm2 results
	assembly$checkm2 <- list(read_checkm2_quality_tsv(fs::path(dir,'checkm2_quality_report.tsv')))
	assembly <- unnest(assembly,checkm2,keep_empty=TRUE,names_sep=".")

	# Add long/short reads stats
	list(
		assembly = assembly,
		contigs = contigs,
		long_reads = samtools_stats_summary(fs::path(dir,'long_reads.bam.stats')),
		short_reads = samtools_stats_summary(fs::path(dir,'short_reads.bam.stats'))
	)
}





#' Read the QC statistics of an assembly written by write_assembly_qc()
#'
#' @param file The JSON file.
#' @return The list of read_assembly_stats(), contigs and depth_hist as tibbles (with all their columns,
#'  even when empty or when all the values are missing).
read_assembly_qc <- function(file) {
	qc <- jsonlite::read_json(file, simplifyVector = TRUE)
	qc$contigs <- as_tibble(qc$contigs) |>
		bind_rows(tibble(
			contig_idx = integer(0), contig_id = character(0), contig_name = character(0), length = numeric(0),
			GC_count = numeric(0), N_count = numeric(0), long_cov = numeric(0), short_cov = numeric(0)
		)) |>
		mutate(across(c(long_cov, short_cov), as.numeric))
	for (r in c("long_reads","short_reads")) {
		if (!is.null(qc[[r]])) {
			qc[[r]]$depth_hist <- as_tibble(qc[[r]]$depth_hist) |>
				bind_rows(tibble(range_max = numeric(0), count = numeric(0)))
		}
	}
	qc
}



#' Table of the QC statistics of assemblies, one row per assembly
#'
#' @param files The JSON files written by write_assembly_qc().
#' @return A tibble with sample_id, assembler_name and the columns assembly.*, long_reads.*, short_reads.*
#'  (without depth_hist) and checkm2.* (NA when missing), sorted by sample_id and assembler_name.
aggregate_assembly_stats <- function(files) {
	fields <- list(
		assembly = c("num_contig","total_len","min_len","median_len","max_len","N50","GC_pct"),
		long_reads = c("num_sequenced","num_mapped","num_mapped_bp","pct_mapped","mean_depth"),
		short_reads = c("num_sequenced","num_mapped","num_mapped_bp","pct_mapped","mean_depth"),
		checkm2 = c("completeness","contamination")
	)
	files |>
		map(function(f) {
			qc <- jsonlite::read_json(f)
			values <- imap(fields, function(names, block) {
				map(names, ~ as.numeric(qc[[block]][[.x]] %||% NA_real_)) |> set_names(paste0(block, ".", names))
			})
			as_tibble(c(list(sample_id = qc$sample_id, assembler_name = qc$assembler_name), list_flatten(unname(values))))
		}) |>
		bind_rows() |>
		arrange(sample_id, assembler_name)
}
