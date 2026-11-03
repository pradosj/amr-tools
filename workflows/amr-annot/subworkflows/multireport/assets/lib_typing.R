
library(tidyverse)
library(GenomicRanges)
library(Biostrings)

read_org_name <- function(txt_file) {
	if (fs::file_exists(txt_file)) {
		readLines(txt_file)	
	} else {
		NA_character_
	}
}


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




read_cgemlst_json <- function(json_file) {
	if (!file.exists(json_file)) return(tibble(mlst_type="?"))
	json <- jsonlite::fromJSON(json_file,simplifyVector=FALSE)
	json$mlst$results$sequence_type |>
		enframe(value = "mlst_type",name = NULL) |>
		mutate(mlst_type=if_else(mlst_type %in% c("Unknown",NA),"?",mlst_type))
}



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



db_load <- function(amr_dir) {

	assemblies <- fs::path(amr_dir,"amr/assemblies") |>
		fs::dir_ls(glob = "*.fasta",fail=FALSE) |>
		enframe(name = NULL,value = "assembly_path") |>
		mutate(assembly_id = fs::path_file(assembly_path) |> fs::path_ext_remove()) |>
		mutate(sample_id = assembly_id |> str_extract("(.*)__(.*)",1)) |>
		mutate(assembler_name = assembly_id |> str_extract("(.*)__(.*)",2)) |>
		mutate(amr_path = assembly_path |> fs::path_ext_set(".amr")) |>
		mutate(org_name = map_chr(fs::path(amr_path,"org_name.txt"),read_org_name)) |>
		mutate(orgfinder_tax = map(fs::path(amr_path,"orgfinder","tax.tsv"),read_orgfinder_tax)) |>
		mutate(orgfinder_ani = map(fs::path(amr_path,"orgfinder","ani.tsv"),read_orgfinder_tsv)) |>
		mutate(contigs = map(assembly_path,contig_meta)) |>
		mutate(mlst = map(fs::path(amr_path,"mlst","data.json"),read_cgemlst_json)) |>
		mutate(plasmidfinder = map(fs::path(amr_path,"plasmidfinder","data.json"),read_plasmidfinder_json)) |>
		mutate(resfinder = map(fs::path(amr_path,"resfinder","data.json"),read_resfinder_json)) |>
		mutate(amrfinderplus = map(fs::path(amr_path,"amrfinderplus","report.tsv"),read_amrfinderplus_tsv)) |>		
		mutate(mobtyper = map(fs::path(amr_path,"mobtyper.tsv"),read_mobtyper_tsv))
	
	samples <- fs::path(amr_dir,"amr/reads") |>
		fs::dir_ls(glob="*.reads_amr",type="directory",fail=FALSE) |>
		enframe(name = NULL,value = "reads_amr_path") |>
		mutate(sample_id = reads_amr_path |> fs::path_file() |> fs::path_ext_remove()) |>
		full_join(select(assemblies,sample_id),by="sample_id") |>
		mutate(resfinder_long_reads = map(fs::path(reads_amr_path,"resfinder_long","data.json"),read_resfinder_json)) |>
		mutate(resfinder_short_reads = map(fs::path(reads_amr_path,"resfinder_short","data.json"),read_resfinder_json)) |>
		mutate(plasmidfinder_long_reads = map(fs::path(reads_amr_path,"plasmidfinder_long","data.json"),read_plasmidfinder_json)) |>
		mutate(plasmidfinder_short_reads = map(fs::path(reads_amr_path,"plasmidfinder_short","data.json"),read_plasmidfinder_json))
	
	list(samples=samples,assemblies=assemblies)
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






