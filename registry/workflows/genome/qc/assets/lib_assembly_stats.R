# QC statistics of an assembly (MERGE_QC_STATS, written as qc.json), its dashboard (genome_qc_report.qmd)
# and the table of all the assemblies (GENOME_MULTIQC). Sourced after lib_io.R.

library(tidyverse)
library(Biostrings)

#' Statistics of each contig of an assembly
#'
#' @param fa_file FASTA file of the assembly.
#' @return A tibble with contig_idx, contig_id (the name up to the first space), contig_name, length,
#'  GC_count and N_count.
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
#'
#' @param stats_file Output of samtools stats (missing file: no reads).
#' @return A one row tibble with num_sequenced, num_mapped, num_mapped_bp, pct_mapped (% of the reads mapped,
#'  NA without reads) and depth_hist (list column: tibble depth, count, the number of bases of the assembly
#'  at each read depth, empty without reads).
samtools_stats_summary <- function(stats_file) {
	s <- read_samtools_stats(stats_file)
	enframe(s) %>% 
	  pivot_wider() %>% 
	  mutate(
	    num_sequenced = as.numeric(pluck(SN, 1, "raw.total.sequences", .default = NA_real_)),
	    num_mapped = as.numeric(pluck(SN, 1, "reads.mapped", .default = NA_real_)),
	    num_mapped_bp = as.numeric(pluck(SN, 1, "bases.mapped..cigar.", .default = NA_real_))
	  ) %>% 
	  mutate(
	    pct_mapped = 100 * num_mapped / num_sequenced,
	    depth_hist = map(COV_HIST,~as_tibble(.) |> select(depth=range_max,count))
	  ) %>% 
	  select(-SN, -COV_HIST)
}

#' Merge all the QC statistics of an assembly found in dir
#'
#' @param dir Folder holding assembly.fasta and, when available, {long,short}_reads.bam.stats (samtools stats),
#'  {long,short}_reads.coverage.tsv (samtools coverage), checkm2_quality_report.tsv and assembly.json (the
#'  asm.json of the assembler, see read_assembly_topology()).
#' @return A list: assembly (one row tibble: num_contig, total_len, min_len, median_len, max_len, N50, GC_pct,
#'  long_meandepth, short_meandepth (mean depth of the contigs weighted by their length, NA without reads),
#'  num_circular (number of circular contigs, NA when the assembler gives no topology),
#'  checkm2.Completeness, checkm2.Contamination, checkm2.Average_Gene_Length (NA without CheckM2)), contigs
#'  (see fa_extract_contigs_stats(), plus long_meandepth, short_meandepth and topology: circular, linear or NA),
#'  long_reads and short_reads
#'  (see samtools_stats_summary()).
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

	# Add the topology of the contigs given by the assembler (circular or linear, NA when unknown)
	contigs <- contigs |>
	  left_join(read_assembly_topology(fs::path(dir,'assembly.json')),by = "contig_id")

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
	    short_meandepth = weighted.mean(short_meandepth,length),
	    num_circular = if (all(is.na(topology))) NA_integer_ else sum(topology == "circular", na.rm = TRUE)
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

#' Read the QC statistics of an assembly (qc.json of MERGE_QC_STATS)
#'
#' @param file The JSON file: the list of merge_genome_qc_stats() and meta (id, assembler_name).
#' @return The same list: meta as a list, assembly, contigs, long_reads and short_reads as tibbles (depth_hist
#'  as a list column of tibbles depth, count), the depth columns numeric and topology
#'  character even when all their values are missing.
read_assembly_qc <- function(file) {
	qc <- jsonlite::read_json(file, simplifyVector = TRUE)
	depth_cols <- c("long_meandepth","short_meandepth")
	qc$assembly <- as_tibble(qc$assembly) |> mutate(across(any_of(depth_cols), as.numeric))
	qc$contigs <- as_tibble(qc$contigs) |>
		bind_rows(tibble(
			contig_idx = integer(0), contig_id = character(0), contig_name = character(0), length = numeric(0),
			GC_count = numeric(0), N_count = numeric(0), long_meandepth = numeric(0), short_meandepth = numeric(0),
			topology = character(0)
		)) |>
		mutate(across(all_of(depth_cols), as.numeric), topology = as.character(topology))
	for (r in c("long_reads","short_reads")) {
		qc[[r]] <- as_tibble(qc[[r]]) |>
			mutate(depth_hist = map(depth_hist, ~ bind_rows(as_tibble(.x), tibble(depth = numeric(0), count = numeric(0)))))
	}
	qc
}

#' Table of the QC statistics of assemblies, one row per assembly
#'
#' @param files The qc.json files of MERGE_QC_STATS.
#' @return A tibble with id, assembler_name, the columns of assembly (checkm2.* included), then the columns of
#'  long_reads and short_reads (without depth_hist) prefixed by long_reads. and short_reads., sorted by id
#'  and assembler_name.
aggregate_assembly_stats <- function(files) {
	files |>
		map(function(f) {
			qc <- read_assembly_qc(f)
			reads <- function(r) qc[[r]] |> select(-depth_hist) |> rename_with(~ paste0(r, ".", .x))
			tibble(id = qc$meta$id, assembler_name = qc$meta$assembler_name) |>
				bind_cols(qc$assembly, reads("long_reads"), reads("short_reads"))
		}) |>
		bind_rows() |>
		arrange(id, assembler_name)
}
