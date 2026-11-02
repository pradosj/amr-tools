
library(GenomicRanges)
library(tidyverse)

# Convert CDS lines of a NCBI GTF into a GFF format suitable to work with bcftools csq
gtf2gff <- function(gtf) {
	cds <- gtf[gtf$type=="CDS"]
	mcols(cds) <- mcols(cds)[c("type","gene_id","gene")]
	cds$Parent <- str_c("transcript:",cds$gene_id)
	cds$phase <- 0
	
	tx <- cds
	tx$ID <- str_c("transcript:",cds$gene_id)
	tx <- stack(range(splitAsList(tx,tx$ID)),"ID")
	tx$type <- "transcript"
	tx$biotype <- "protein_coding"
	tx$Parent <- str_replace(as.character(tx$ID),"^transcript:","gene:")
	
	genes <- tx
	genes$ID <- str_replace(as.character(tx$ID),"^transcript:","gene:")
	genes$type <- "gene"
	genes$biotype <- "protein_coding"
	genes$Parent <- NULL
	#genes$Name <- genes$gene
	
	gene_names <- mcols(gtf) %>% 
		as_tibble() %>% 
		mutate(Name = gene) %>% 
		filter(!is.na(Name)) %>% 
		mutate(ID=str_c("gene:",gene_id)) %>% 
		select(ID,Name) %>% distinct() %>% 
		group_by(ID) %>% 
		slice_head(n=1) %>% 
		ungroup()

	gff <- c(genes,tx,cds)
	gff$gene_id <- gff$gene <- NULL
	gff$Name <- gene_names$Name[match(gff$ID,gene_names$ID)]
	gff
}


