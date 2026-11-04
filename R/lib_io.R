# Readers of the files of the bioinformatics tools and of the common formats, shared by the R steps of the
# pipelines (reached through file links). It attaches no package, so that the order of the packages set
# by each caller is kept: the callers attach tidyverse, Biostrings (FASTA files) and tidygraph (read_tax()).


# ---------------------------------------
#  Common formats
# ---------------------------------------

#' Read a BLAST tabular output
#'
#' @param f BLAST output, as produced with `-outfmt "6 std qlen slen"`.
#' @return A tibble with one row per alignment.
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

#' Read a FASTA file
#'
#' @param f FASTA file.
#' @return A DNAStringSet, its sequences named by their id (first word of the header), as BLAST does.
read_fasta <- function(f) {
	fa <- readDNAStringSet(f)
	names(fa) <- sub("\\s.*", "", names(fa))
	fa
}

#' Contigs of a FASTA file, with the tags of their header
#'
#' The tags are the `[key=value]` or `key=value` modifiers of the headers, e.g. `[topology=circular]`.
#'
#' @param fasta_filename FASTA file.
#' @return A tibble with one row per contig: contig_id, contig_length, GC, title, and one tag_<key> column per tag.
contig_meta <- function(fasta_filename) {
	fasta_filename <- as.character(fasta_filename)
	fa <- Biostrings::readDNAStringSet(fasta_filename)
	meta <- tibble(
		contig_id = str_replace(names(fa)," .*",""),
		contig_length = lengths(fa),
		header = names(fa),
		GC = as.vector(letterFrequency(fa,"GC",as.prob=TRUE))
	)
	stopifnot("FASTA contains duplicated contig_id" = all(!duplicated(meta$contig_id)))
	
	# regex-pattern for modifier and header line
	tag_pat <- "((\\[([^=]+)=([^\\]]+)\\])|(([^=]+)=([^\\s]+)))"
	hdr_pat <- str_glue("^([^\\s]+)((\\s*{tag_pat})*)\\s*(.*)$")
	if (!all(str_detect(meta$header,hdr_pat))) rlang::abort("Some header lines are malformed")
	
	meta <- meta |>
		mutate(title = str_extract(header,hdr_pat,group=11)) |>
		mutate(tags_str=str_trim(str_extract(header,hdr_pat,group=2)))
	
	
	# Extract key-value pairs from modifiers
	tags <- tibble(select(meta,contig_id,tags_str)) |>
		mutate(tags = str_extract_all(tags_str,tag_pat),tags_str = NULL) |>
		unnest(tags) |>
		mutate(tags = str_trim(tags)) |>
		# TODO: use coalesce() here ?
		mutate(
			tag_name = if_else(is.na(str_extract(tags,tag_pat,3)),str_extract(tags,tag_pat,6),str_extract(tags,tag_pat,3)),
			tag_value = if_else(is.na(str_extract(tags,tag_pat,4)),str_extract(tags,tag_pat,7),str_extract(tags,tag_pat,4)),
			tags = NULL
		)
	
	stopifnot("Some modifiers are set several times in the same header" = !(summarize(tags,.by=c(contig_id,tag_name),any(n()>1)) |> pull() |> any()))
	expected_tags_structure <- tibble(
		tag_topology = character(0)
	)
	tags <- pivot_wider(tags,id_cols="contig_id",names_from = "tag_name",values_from = "tag_value",names_prefix = "tag_") |>
		bind_rows(expected_tags_structure)
	
	left_join(meta,tags,by="contig_id",relationship = "one-to-one") |>
		select(!c(tags_str,header))
}


# ---------------------------------------
#  AMR annotation tools (AMR_ANNOT)
# ---------------------------------------

#' Organism name found by orgfinder
#'
#' @param txt_file org_name.txt file of an assembly.
#' @return The organism name, or NA when the file does not exist.
read_org_name <- function(txt_file) {
	if (fs::file_exists(txt_file)) {
		readLines(txt_file)	
	} else {
		NA_character_
	}
}

#' Hits of a ResFinder result
#'
#' @param json_file ResFinder data.json.
#' @return A tibble with one row per hit (empty when the file does not exist).
read_resfinder_json <- function(json_file) {
	expected_structure <- tibble(
		contig_id = character(0),
		name = character(0),
		query_start_pos = numeric(0),
		query_end_pos = numeric(0),
		query_id = character(0),
		coverage = numeric(0),
		identity = numeric(0)
	)
	if (fs::file_exists(json_file)) {
		out <- jsonlite::fromJSON(file(json_file,"rb"),simplifyVector=FALSE) |>
			pluck("seq_regions") |>
			tibble() |> 
			unnest_wider(1) |>
			bind_rows(tibble(query_id = character(0))) |>
			mutate(contig_id = str_replace(query_id," .*",""))
	} else {
		out <- NULL
	}
	out |>
		bind_rows(expected_structure)	|>
		mutate(position=str_glue("{query_start_pos}-{query_end_pos}")) |>
		relocate(contig_id,resistance_name=name,query_start_pos,query_end_pos,coverage,identity,position)
}

#' Hits of an AMRFinderPlus result
#'
#' @param tsv_file AMRFinderPlus report.tsv.
#' @return A tibble with one row per hit (empty when the file does not exist).
read_amrfinderplus_tsv <- function(tsv_file) {
	expected_structure <- tibble(
		contig_id = character(0),
		resistance_name = character(0),
		coverage = numeric(0),
		identity = numeric(0),
		position = character(0)
	)	
	if (!file.exists(tsv_file)) return(expected_structure)
	read_tsv(tsv_file,show_col_types = FALSE,col_types = cols(`Contig id`="c",`% Coverage of reference`="n",`% Identity to reference`="n",`Element symbol`="c",.default="c")) |>
		dplyr::rename(contig_id=`Contig id`,resistance_name=`Element symbol`,coverage=`% Coverage of reference`,identity=`% Identity to reference`) |>
		mutate(position=str_glue("{Start}-{Stop}:{Strand}")) |>
		bind_rows(expected_structure) |>
		relocate(contig_id,coverage,identity,resistance_name,position)
}

#' Hits of a PlasmidFinder result
#'
#' @param json_file PlasmidFinder data.json.
#' @return A tibble with one row per hit (empty when the file does not exist).
read_plasmidfinder_json <- function(json_file) {
	expected_structure <- tibble(
		contig_name = character(0),
		plasmid = character(0),
		coverage = numeric(0),
		identity = numeric(0),
		positions_in_contig = character(0)
	)
	if (fs::file_exists(json_file)) {
		out <- jsonlite::fromJSON(file(json_file,"rb"),simplifyVector=FALSE) |>
			pluck("plasmidfinder","results") |>
			enframe(name = "db_lev1") |>
			unnest_longer(value,indices_to = "db_lev2") |>
			filter(value != "No hit found") |>		
			unnest_longer(value,indices_include = FALSE) |>
			unnest_wider(value)
	} else {
		out <- NULL
	}
	out |>
		bind_rows(expected_structure) |>
		mutate(contig_id = str_replace(contig_name," .*","")) |>
		relocate(contig_id,plasmid_type=plasmid,coverage,identity)
}

#' Sequence type of a CGE MLST result
#'
#' @param json_file CGE MLST data.json.
#' @return A tibble with the column mlst_type ("?" when unknown or when the file does not exist).
read_cgemlst_json <- function(json_file) {
	if (!file.exists(json_file)) return(tibble(mlst_type="?"))
	json <- jsonlite::fromJSON(json_file,simplifyVector=FALSE)
	json$mlst$results$sequence_type |>
		enframe(value = "mlst_type",name = NULL) |>
		mutate(mlst_type=if_else(mlst_type %in% c("Unknown",NA),"?",mlst_type))
}

#' Replicon and relaxase types of a MOB-typer result
#'
#' @param tsv_file MOB-typer mobtyper.tsv.
#' @return A tibble with one row per contig (empty when the file does not exist).
read_mobtyper_tsv <- function(tsv_file) {
		expected_structure <- tibble(
			contig_id = character(0),
			contig_name = character(0),
			relaxase_types = character(0),
			rep_types = character(0)
		)	
		if (!file.exists(tsv_file)) return(expected_structure)
		read_tsv(tsv_file,col_types = cols(.default="c")) |>
			mutate(contig_id=str_replace(sample_id," .*","")) |>
			dplyr::rename(contig_name=sample_id,relaxase_types=`relaxase_type(s)`,rep_types=`rep_type(s)`) |>
			bind_rows(expected_structure) |>
			relocate(contig_id)
}

#' Taxonomy of the organisms of an orgfinder result
#'
#' @param tsv_file orgfinder tax.tsv.
#' @return A tibble with one row per organism (empty when the file does not exist).
read_orgfinder_tax <- function(tsv_file) {

	expected_structure <- tibble(
		assembly_acc = character(0),
		org_name = character(0),
		tax_id = character(0),
		species_id = character(0),
		genus_id = character(0),
		family_id = character(0),
		order_id = character(0),
		class_id = character(0),
		phylum_id = character(0),
		species_name = character(0),
		genus_name = character(0),
		family_name = character(0),
		order_name = character(0),
		class_name = character(0),
		phylum_name = character(0)
	)
	if (!file.exists(tsv_file)) return(expected_structure)
	read_tsv(tsv_file,show_col_types = FALSE,col_types = cols(.default = "c")) |>
		bind_rows(expected_structure) |>
		mutate(species_name = if_else(is.na(species_name),org_name,species_name))
}

#' FastANI hits of an orgfinder result
#'
#' @param tsv_file orgfinder ani.tsv.
#' @return A tibble with one row per reference genome hit (empty when the file does not exist).
read_orgfinder_tsv <- function(tsv_file) {
	expected_structure <- tibble(
		query = character(0),
		assembly_acc = character(0),
		ANI = numeric(0),
		bi_frag = integer(0),
		query_frag = integer(0),		
		org_name = character(0),
		tax_id = character(0)
	)	
	if (!file.exists(tsv_file)) return(expected_structure)
	read_tsv(tsv_file,show_col_types = FALSE,col_types = cols(ANI = "n",bi_frag="i",query_frag="i",.default = "c")) |>
		bind_rows(expected_structure)
}



# ---------------------------------------
#  Assembly QC tools (ASSEMBLIES_QC)
# ---------------------------------------

#' samtools stats output
#'
#' @param f samtools stats file.
#' @return A list with SN (the summary numbers, as a named list) and COV_HIST (the coverage histogram);
#'   empty when the file does not exist.
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

#' Completeness and contamination of a CheckM2 result
#'
#' @param tsv_file CheckM2 quality_report.tsv.
#' @return A tibble with Completeness, Contamination and Average_Gene_Length (empty when the file does not exist).
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


# ---------------------------------------
#  NCBI taxonomy (orgfinder database)
# ---------------------------------------

#' NCBI taxonomy
#'
#' @param taxdump_dir Folder of the NCBI taxdump (names.dmp, nodes.dmp).
#' @return A tbl_graph of the taxa (tax_id, parent, rank, tax_name), from the parents to their children.
read_tax <- function(taxdump_dir) {
	SN <- read_tsv(fs::path(taxdump_dir,"names.dmp"),col_names = c("tax_id","tax_name","name_class"),col_types="c_c___c_") |> 
		filter(name_class=="scientific name") |>
		mutate(name_class=NULL)
	N <- read_tsv(fs::path(taxdump_dir,"nodes.dmp"),col_names = c("tax_id","parent","rank"),col_types="c_c_c_____________________") |>
		left_join(SN,by="tax_id",relationship="one-to-one")
	tax <- tbl_graph(N,select(N,parent,tax_id),node_key = "tax_id")
	tax
}
