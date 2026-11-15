

library(tidyverse)
library(Biostrings)

fa_extract_contigs_stats <- function(fa_file) {
	fa <- readDNAStringSet(fa_file)
	contigs <- tibble(
		contig_idx = seq_along(fa),
		contig_name = str_replace(names(fa)," .*",""),
		length = lengths(fa),
		GC_count = as.vector(Biostrings::letterFrequency(fa,"GC")),
		N_count = as.vector(Biostrings::letterFrequency(fa,"N"))
	)
	contigs
}


read_assembly_stats <- function(dir) {
	stats <- list()
	stats$contigs <- fa_extract_contigs_stats(fs::path(dir,'assembly.fasta'))
	stats$short_reads <- read_samtools_stats(fs::path(dir,'short_reads.bam.stats'))
	stats$long_reads <- read_samtools_stats(fs::path(dir,'long_reads.bam.stats'))
	stats$checkm2 <- read_checkm2_quality_tsv(fs::path(dir,'checkm2_quality_report.tsv'))
	
	stats$assembly_stats <- stats$contigs |>
		summarize(
			num_contig = n(),
			total_len = sum(length),
			min_len = min(length),
			median_len = median(length),
			max_len = max(length),
			N50 = N50(length),
			GC_pct = sum(GC_count) / total_len
		)
	stats
}

aggregate_assembly_stats <- function(...) {
	x <- list(...) |>
		enframe(name = "assembly_idx") |>
		unnest_wider(value) 
	if (nrow(x)<1) return(x)
	x |>
		select(!contigs) |>
		mutate(
			across(c(short_reads,long_reads),~tibble(
				num_sequenced = map_dbl(.x,~pluck(.x,"SN","raw.total.sequences",.default=NA_real_)),
				num_mapped = map_dbl(.x,~pluck(.x,"SN","reads.mapped",.default=NA_real_)),
				num_mapped_bp = map_dbl(.x,~pluck(.x,"SN","bases.mapped..cigar.",.default=NA_real_)),
				pct_mapped = num_mapped/num_sequenced,
				avg_coverage_depth = num_mapped_bp / assembly_stats$total_len
			))
		) |>
		unnest_wider(c(assembly_stats,short_reads,long_reads,checkm2),names_sep = ".") |>
		relocate("assembly_idx","meta",starts_with("assembly_stats"))
}

