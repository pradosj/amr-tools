#!/usr/bin/env Rscript

db_dir <- "/Users/prados/Documents/AMR-genomics/amr-tools/amr-tools/pipelines/orgfinder/results/db/orgfinder"
fastani_file <- "/Users/prados/Documents/AMR-genomics/amr-tools/amr-tools/pipelines/orgfinder/results/GCF_011397115.1.fastani"
	
db <- readr::read_tsv(fs::path(db_dir,"db.tsv"),col_types = "ccc")

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







