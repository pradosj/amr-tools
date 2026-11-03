#!/usr/bin/env Rscript

db_tsv <- "./db/db.tsv"
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
readr::write_tsv(out,file='orgfinder/ani.tsv')

# Best hit (highest ANI): accession, organism name and ANI, as shell variables
best <- dplyr::slice_max(out,ANI,n=1,with_ties=FALSE)
writeLines(c(
	sprintf('ORG_ACC="%s"',c(best$assembly_acc,"")[1]),
	sprintf('ORG_NAME="%s"',c(best$org_name,"")[1]),
	sprintf('ORG_ANI="%s"',c(best$ANI,"")[1])
),"org.env")
