#!/usr/bin/env Rscript

base_dir <- './'
ncbi_dir <- fs::path(base_dir,"ncbi_db")
db_dir <- fs::path(base_dir,"orgfinder_db")

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Load the list of reference genome downloaded from NCBI
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
library(tidyverse)
ss <- read_tsv(file.path(ncbi_dir,"accessions.tsv"),col_types = cols("Organism Taxonomic ID"="c")) |>
	left_join(
		list.files(ncbi_dir,"_genomic.fna$",recursive = TRUE,full.names = TRUE) |>
			enframe(value = "src_path",name=NULL) |>
			mutate(`Assembly Accession`=basename(dirname(src_path)))
	) |>
	mutate(db_path=str_glue('fna/{`Assembly Accession`}.fna')) |>
	group_by(`Assembly Accession`) |>
	slice_head(n=1) |>
	ungroup()



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# TODO: Add molecule type and exclude plasmids
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# function() {
# 	mol_type <- read_tsv(file.path(ncbi_dir,"molecule_types.tsv"),col_names = c("Assembly Accession","contig","type"),col_types = "c") %>%
# 		filter(type != "Plasmid")
# 	ss <- nest_join(ss,mol_type,by="Assembly Accession")
# 	ss %>%
# 		rowwise() %>%
# 		mutate({
# 			dna <- Biostrings::readDNAStringSet(src_path)
# 			names(dna) <- str_replace(names(dna)," .*","")
# 			all(names(dna) %in% mol_type$contig)
# 		})
# }
	


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Copy the FASTA to the DB folder with appropriate renaming
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
dir.create(file.path(db_dir,"fna"),recursive = TRUE)
file.copy(ss$src_path,file.path(db_dir,ss$db_path),overwrite = TRUE)

# Build db annotation table
message("generate db.tsv")
db <- ss |>
	select(assembly_acc=`Assembly Accession`,org_name=`Organism Name`,tax_id=`Organism Taxonomic ID`)
write_tsv(db,file = file.path(db_dir,"db.tsv"),na="")



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Prepare Taxonomy
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
library(tidygraph)

# Load NCBI taxonomy
read_tax <- function(taxdump_dir) {
	SN <- read_tsv(fs::path(taxdump_dir,"names.dmp"),col_names = c("tax_id","tax_name","name_class"),col_types="c_c___c_") |> 
		filter(name_class=="scientific name") |>
		mutate(name_class=NULL)
	N <- read_tsv(fs::path(taxdump_dir,"nodes.dmp"),col_names = c("tax_id","parent","rank"),col_types="c_c_c_____________________") |>
		left_join(SN,by="tax_id",relationship="one-to-one")
	tax <- tbl_graph(N,select(N,parent,tax_id),node_key = "tax_id")
	tax
}

message("load taxonomy")
tax <- read_tax(fs::path(base_dir,"taxdump")) |>
	mutate(is_db_tax = tax_id %in% db$tax_id) 

# Find all ancestors of selected nodes
message("find ancestors")
db_ancestors <- igraph::ego(tax,order=igraph::vcount(tax),nodes=igraph::V(tax)[is_db_tax],mode="in") |> 
	map(~.x$tax_id) |>
	setNames(igraph::V(tax)[is_db_tax]$tax_id) |>
	enframe(name = "leaf_id",value = "ancestor_id") |>
	unnest(ancestor_id)

# Subset the taxonomy to selected elements and its ancestors
message("subset taxonomy to selected elements")
TAX <- tax |>
	activate(edges) |>
	filter(!edge_is_loop()) |>
	activate(nodes) |>
	mutate(is_db_ancestor = tax_id %in% db_ancestors$ancestor_id) |>
	filter(is_db_ancestor)
saveRDS(TAX,file.path(db_dir,"tax.rds"))

# Add genus informations to the database
DB <- db |>
	left_join(
		inner_join(
			select(db_ancestors,tax_id=leaf_id,id=ancestor_id),
			select(as_tibble(TAX),id=tax_id,rank=rank,name=tax_name) |>
				filter(rank %in% c("superkingdom","phylum","class","order","family","genus","species"))
		) |>
			pivot_wider(id_cols = "tax_id",names_from = "rank",values_from = c("id","name"),names_glue = "{rank}_{.value}")
	)
write_tsv(DB,file.path(db_dir,"tax.tsv"),na = "")





