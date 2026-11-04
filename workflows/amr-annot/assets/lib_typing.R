
library(tidyverse)
library(GenomicRanges)
library(Biostrings)

#' AMR results of assemblies
#'
#' @param amr_path Results folders of the assemblies (amr/assemblies/<assembly_id>.amr).
#' @param assembly_path FASTA files of the assemblies.
#' @param sample_id,assembler_name Sample id and assembler name of the assemblies.
#' @return A tibble with one row per assembly, the results of each tool in a list column.
read_assemblies_amr <- function(amr_path, assembly_path, sample_id, assembler_name) {
	tibble(amr_path = fs::path(amr_path), assembly_path = fs::path(assembly_path), sample_id, assembler_name) |>
		mutate(assembly_id = str_c(sample_id,"__",assembler_name)) |>
		mutate(org_name = map_chr(fs::path(amr_path,"org_name.txt"),read_org_name)) |>
		mutate(orgfinder_tax = map(fs::path(amr_path,"orgfinder","tax.tsv"),read_orgfinder_tax)) |>
		mutate(orgfinder_ani = map(fs::path(amr_path,"orgfinder","ani.tsv"),read_orgfinder_tsv)) |>
		mutate(contigs = map(assembly_path,contig_meta)) |>
		mutate(mlst = map(fs::path(amr_path,"mlst","data.json"),read_cgemlst_json)) |>
		mutate(plasmidfinder = map(fs::path(amr_path,"plasmidfinder","data.json"),read_plasmidfinder_json)) |>
		mutate(resfinder = map(fs::path(amr_path,"resfinder","data.json"),read_resfinder_json)) |>
		mutate(amrfinderplus = map(fs::path(amr_path,"amrfinderplus","report.tsv"),read_amrfinderplus_tsv)) |>
		mutate(mobtyper = map(fs::path(amr_path,"mobtyper.tsv"),read_mobtyper_tsv)) |>
		select(!c(amr_path,assembly_path)) |>
		relocate(assembly_id,sample_id,assembler_name)
}

#' AMR results of reads samples
#'
#' @param reads_amr_path Results folders of the samples (amr/reads/<sample_id>.reads_amr).
#' @param sample_id Sample ids.
#' @return A tibble with one row per sample, the results of each tool in a list column.
read_reads_amr <- function(reads_amr_path, sample_id) {
	tibble(reads_amr_path = fs::path(reads_amr_path), sample_id) |>
		mutate(resfinder_long_reads = map(fs::path(reads_amr_path,"resfinder_long","data.json"),read_resfinder_json)) |>
		mutate(resfinder_short_reads = map(fs::path(reads_amr_path,"resfinder_short","data.json"),read_resfinder_json)) |>
		mutate(plasmidfinder_long_reads = map(fs::path(reads_amr_path,"plasmidfinder_long","data.json"),read_plasmidfinder_json)) |>
		mutate(plasmidfinder_short_reads = map(fs::path(reads_amr_path,"plasmidfinder_short","data.json"),read_plasmidfinder_json)) |>
		select(!reads_amr_path)
}

#' AMR results of all the assemblies and reads samples
#'
#' @param assemblies_rds RDS files saved from read_assemblies_amr().
#' @param reads_rds RDS files saved from read_reads_amr().
#' @return A list with the tibbles `samples` and `assemblies`.
db_aggregate <- function(assemblies_rds, reads_rds) {
	list(
		samples = bind_rows(read_reads_amr(character(0),character(0)),map(reads_rds,readRDS)),
		assemblies = bind_rows(read_assemblies_amr(character(0),character(0),character(0),character(0)),map(assemblies_rds,readRDS))
	)
}

#' AMR results of a published amr/ folder
#'
#' @param amr_dir The folder holding amr/.
#' @return A list with the tibbles `samples` and `assemblies`, as db_aggregate().
db_load <- function(amr_dir) {
	assembly_path <- fs::dir_ls(fs::path(amr_dir,"amr/assemblies"),glob = "*.fasta",fail=FALSE)
	assembly_id <- fs::path_file(assembly_path) |> fs::path_ext_remove()
	reads_amr_path <- fs::dir_ls(fs::path(amr_dir,"amr/reads"),glob="*.reads_amr",type="directory",fail=FALSE)
	list(
		samples = read_reads_amr(reads_amr_path,fs::path_file(reads_amr_path) |> fs::path_ext_remove()),
		assemblies = read_assemblies_amr(
			fs::path_ext_set(assembly_path,".amr"),
			assembly_path,
			str_extract(assembly_id,"(.*)__(.*)",1),
			str_extract(assembly_id,"(.*)__(.*)",2)
		)
	)
}

#' Summary tables of the AMR results
#'
#' @param db AMR results, as returned by db_aggregate() or db_load().
#' @return A list of tibbles: assemblies, contigs, resistances, plasmidfinder_hits.
amr_tables <- function(db) {
	tbl <- list(
		assemblies = summarise_assembly(db),
		contigs = summarise_contigs(db),
		resistances = summarise_resistances(db),
		plasmidfinder_hits = summarise_plasmidfinder_hits(db)
	)
	tbl$contigs <- left_join(tbl$contigs,select(tbl$assemblies,sample_id,assembler_name,species_name,mlst_type),by=c('sample_id','assembler_name'))
	tbl
}


summarise_assembly <- function(db) {
	
	assemblies <- db$assemblies |> 
		select(assembly_id,sample_id,assembler_name,contigs) |> 
		unnest(contigs) |> 
		group_by(assembly_id,sample_id,assembler_name) |> 
		bind_rows(tibble(contig_length=numeric(0),GC=numeric(0))) |>
		summarise(num_contig=n(),assembly_length=sum(contig_length),GC=weighted.mean(GC,contig_length),N50=N50(contig_length)) |>
		ungroup()
	mlst <- db$assemblies |> select(assembly_id,mlst) |> unnest(mlst) |> bind_rows(tibble(mlst_type=character(0))) |> select(assembly_id,mlst_type)

	org_name <- db$assemblies |> 
		select(assembly_id,org_name,orgfinder_tax) |>
		mutate(org=map2(orgfinder_tax,org_name,~left_join(enframe(.y,name = NULL,value="org_name"),.x) |> select(!org_name) )) |>
		unnest(org) |>
		bind_rows(tibble(species_name=character(0),genus_name=character(0))) |>
		select(assembly_id,org_name,species_name,genus_name)

	orgfinder <- db$assemblies |> 
		mutate(orgfinder = map2(orgfinder_ani,orgfinder_tax,~left_join(.x,.y,by="org_name"))) |>
		select(assembly_id,orgfinder) |> 	
		unnest(orgfinder) |> 
		bind_rows(tibble(org_name=character(0),ANI=numeric(0),species_name=character(0),genus_name=character(0))) |>
		select(assembly_id,org_name,ANI,species_name,genus_name) |>
		group_by(assembly_id) |>
		slice_max(ANI,n=1,with_ties = FALSE) |>
		ungroup()
	
	assemblies  |>
		relocate(assembly_id) |>
		left_join(org_name,by="assembly_id",relationship = "one-to-one")	|>
		left_join(mlst,by="assembly_id",relationship = "one-to-one") |>
		left_join(rename_with(orgfinder,.cols=!assembly_id,~str_c("orgfinder.",.)),by="assembly_id",relationship = "one-to-one")
}

summarise_resistances <- function(db) {
	bind_rows(
		resfinder = db$assemblies |> 
			select(assembly_id,sample_id,assembler_name,resfinder) |> 
			unnest(resfinder) |>
			bind_rows(tibble(contig_id=character(0),resistance_name=character(0),coverage=numeric(0),identity=numeric(0),position=character(0))) |>
			select(assembly_id,sample_id,assembler_name,contig_id,resistance_name,coverage,identity,position),
		amrfinderplus = db$assemblies |> 
			select(assembly_id,sample_id,assembler_name,amrfinderplus) |> 
			unnest(amrfinderplus) |>
			bind_rows(tibble(contig_id=character(0),resistance_name=character(0),coverage=numeric(0),identity=numeric(0),position=character(0))) |>
			select(assembly_id,sample_id,assembler_name,contig_id,resistance_name,coverage,identity,position),
		resfinder_long_reads = db$samples |> 
			select(sample_id,resfinder_long_reads) |> 
			unnest(resfinder_long_reads) |>
			bind_rows(tibble(contig_id=character(0),resistance_name=character(0),coverage=numeric(0),identity=numeric(0),position=character(0))) |>
			select(sample_id,contig_id,resistance_name,coverage,identity,position),
		resfinder_short_reads = db$samples |> 
			select(sample_id,resfinder_short_reads) |> 
			unnest(resfinder_short_reads) |>
			bind_rows(tibble(contig_id=character(0),resistance_name=character(0),coverage=numeric(0),identity=numeric(0),position=character(0))) |>
			select(sample_id,contig_id,resistance_name,coverage,identity,position),
		.id = "source"
	) |>
		arrange(assembly_id,sample_id,assembler_name,contig_id,resistance_name,desc(coverage),desc(identity)) |>
		ungroup()
}

summarise_plasmidfinder_hits <- function(db) {
	select(db$assemblies,assembly_id,sample_id,assembler_name,plasmidfinder) |> 
		unnest(plasmidfinder) |>
		bind_rows(tibble(contig_id=character(0),plasmid_type=character(0),coverage=numeric(0),identity=numeric(0),positions_in_contig=character(0))) |>
		select(assembly_id,sample_id,assembler_name,contig_id,plasmid_type,coverage,identity,position=positions_in_contig) |>
		left_join(summarise_contigs(db) |> select(assembly_id,contig_id,contig_length,contig_is_plasmid=is_plasmid)) |>
		arrange(assembly_id,sample_id,assembler_name,contig_id,desc(coverage),desc(identity)) |>
		ungroup()
}

summarise_contigs <- function(db) {
	contigs <- db$assemblies |> 
		select(assembly_id,sample_id,assembler_name,contigs) |> 
		unnest(contigs)
	res <- summarise_resistances(db) |>
		group_by(assembly_id,contig_id) |>
		summarise(resistance_names = list(str_unique(resistance_name)))
	plf <- db$assemblies |> 
		select(assembly_id,plasmidfinder) |> 
		unnest(plasmidfinder) |>
		bind_rows(tibble(contig_id=character(0),plasmid_type=character(0))) |>
		group_by(assembly_id,contig_id) |>
		summarise(plasmid_types=list(str_unique(plasmid_type)))
	mob <- db$assemblies |> 
		select(assembly_id,mobtyper) |> 
		unnest(mobtyper) |>
		bind_rows(tibble(contig_id=character(0),relaxase_types=character(0))) |>
		mutate(relaxase_types = str_split(relaxase_types,",") |> map(str_unique)) |>
		select(assembly_id,contig_id,relaxase_types)
	contigs |>
		bind_rows(tibble(contig_id=character(0),contig_length=numeric(0),GC=numeric(0),tag_topology=character(0))) |>
		select(assembly_id,sample_id,assembler_name,contig_id,contig_length,GC,topology=tag_topology) |>
		left_join(res,by=c("assembly_id","contig_id"),relationship = "one-to-one") |>
		left_join(plf,by=c("assembly_id","contig_id"),relationship = "one-to-one") |>
		left_join(mob,by=c("assembly_id","contig_id"),relationship = "one-to-one")  |>
		mutate(is_plasmid = (contig_length<=500000) & lengths(plasmid_types)>0)
}






