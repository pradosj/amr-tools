

library(tidyverse)
library(Biostrings)

read_samtools_stats <- function(f) {
	if (fs::file_exists(f)) {
		stats <- readLines(f)	
	} else {
		stats <- character(0)
	}
	list(
		SN = read.table(
				text = str_subset(stats,"^SN\t"),
				sep = "\t",comment="",fill = TRUE,
				col.names = c("tag","var","value","comment")
			) |>
			mutate(var=str_replace(var,": *$","") |> make.names()) |>
			pull(value,name = var) |>
			as.list(),
		COV_HIST = read.table(
					text = str_subset(stats,"^COV\t"),
					sep = "\t",comment="",fill = TRUE,
					col.names = c("tag","range","range_max","count"),colClasses = c(tag="NULL")
				) 
	)
}

read_checkm2_quality_tsv <- function(tsv_file) {
	#tsv_file <- "results/assemblies/RH2_custom_assembly1__none.qc/checkm2_quality.tsv"
	expected_structure <- tibble(
		Completeness = numeric(0),
		Contamination = numeric(0),
		Average_Gene_Length = numeric(0)
	)	
	if (!file.exists(tsv_file)) return(expected_structure)
	read_tsv(tsv_file,show_col_types = FALSE,col_types = cols(Completeness="n",Contamination="n",Average_Gene_Length="n",.default="-")) |>
		bind_rows(expected_structure) |>
		relocate(Completeness,Contamination,Contamination)
}

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

