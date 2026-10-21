#!/usr/bin/env Rscript

db_tsv <- "./db.tsv"
fastani_file <- "./input.fastani"
	
db <- readr::read_tsv(db_tsv,col_types = "ccc")

# Load FastANI output and add DB information
out <- readr::read_tsv(fastani_file,col_names = c("query","ref","ANI","bi_frag","query_frag")) |>
	dplyr::mutate(
		query = fs::path_file(query),
		assembly_acc = fs::path_ext_remove(fs::path_file(ref)),
		ref = NULL
	) |>
	dplyr::left_join(db,by="assembly_acc",relationship="many-to-one") |>
	dplyr::relocate(query,assembly_acc) 
readr::write_tsv(out,file='output.orgfinder')







