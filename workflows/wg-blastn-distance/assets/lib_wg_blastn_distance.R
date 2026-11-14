source("assets/lib_io.R")
source("assets/lib_blast.R")



#' Distance graph of the sequences, from their BLASTN against themselves
#'
#' For each pair of sequences, sum mismatches, gapopen, and extract number of fragments.
#' And keep a single 
#'
#' @param aln Blast alignments, as read by read_blast6().
#' @return An undirected tbl_graph: the sequences as nodes (name), and an edge per pair of sequences aligning (each
#'   sequence with itself included), with num_mis, num_gapopen, n_frag, total_frag_len, err and weight.
blast_undirected_distance_graph <- function(aln) {
  g <- aln |>
  
  tbl_graph(
    nodes = tibble(name = unique(c(g$qid, g$sid))),
    edges = g,
    directed = FALSE
  )
}



#' Distance graph of genomes, from the BLASTN of all their sequences against themselves
#'
#' The alignments of at least min_len bp that are not included in another alignment of the same query on the same
#' subject region (see blast_filter_inner_alignments()) are pooled by pair of genomes (all the sequences of a genome
#' together), then scored with blast_distance_graph(): the weight of an edge is the number of errors per aligned bp,
#' the worst of the two directions. The unaligned parts of the genomes are not counted.
#'
#' @param blast Self BLASTN of the sequences of all the genomes, in the format of read_blast6().
#' @param seq_genomes The genome of each sequence: columns meta_id (the genome) and new_seq_id (the sequence id in the
#'   BLAST), as written by FASTA_CAT (cat.tsv).
#' @param min_len Minimum length of the alignments.
#' @return An undirected tbl_graph (see blast_distance_graph()): the aligning genomes as nodes (name), an edge per pair
#'   of genomes aligning (each genome with itself included).
genome_distance_graph <- function(blast, seq_genomes, min_len = 1000) {

  hits <- read_blast6(blast) |>
		filter(align_len >= min_len) |>
		blast_filter_inner_alignments()
  
  # undirected seq 2 seq distance summary
  seq2seq <- hits |>
    group_by(query_id, subject_id) |>
    summarise(
      num_mis = sum(num_mismatch),
      num_gapopen = sum(num_gapopen),
      n_frag = n(),
      total_frag_len = sum(align_len)
    ) |>
    mutate(
      err = (num_mis + (n_frag - 1)),
      weight = err / total_frag_len
    ) |>
    group_by(
      eid1 = pmin(as.character(query_id), as.character(subject_id)),
      eid2 = pmax(as.character(query_id), as.character(subject_id))
    ) |>
    slice_max(weight, n = 1L, with_ties = FALSE) |>
    ungroup() |>
    left_join(select(seq_genomes, eid1 = new_seq_id, genome1 = meta_id), by="query_id") |> 
    left_join(select(seq_genomes, eid2 = new_seq_id, genome2 = meta_id), by="subject_id") |>
    select(eid1,eid2,genome1,genome2,num_mis,num_gapopen,n_frag,total_frag_len,err,weight)

  # keep the best inter genome reciprocal hit
  seq2seq <- seq2seq |> 
    semi_join(
      by = c("eid1","eid2"),
      bind_rows(
        select(seq2seq,eid1,eid2,seq=eid1,genome=genome2,total_frag_len),
        select(seq2seq,eid1,eid2,seq=eid2,genome=genome1,total_frag_len)
      ) |>
      group_by(seq,genome) |>
      slice_max(total_frag_len,n=1,with_ties = FALSE),
    )
  
  # summarize per genome
  gen2gen <- seq2seq %>% 
    group_by(
      g1 = pmin(as.character(genome1), as.character(genome2)),
      g2 = pmax(as.character(genome1), as.character(genome2))
    ) |>
    summarise(
      num_mis = sum(num_mismatch),
      num_gapopen = sum(num_gapopen),
      n_frag = sum(n_frag),
      total_frag_len = sum(total_frag_len),
      err = sum(err)
    ) |>
    mutate(
      weight = err / total_frag_len
    )
  
  tbl_graph(
    nodes = tibble(name = unique(seq_genomes$meta_id)),
    edges = gen2gen,
    directed = FALSE
  )
}
