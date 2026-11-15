# wg-blastn-distance

`wg-blastn-distance` a nextflow pipeline computing the distances between whole genomes, without annotation
(`WG_BLASTN_DISTANCE`, see [`wg-blastn-distance`](../../../workflows/evo/wg-blastn-distance)):

1. all the sequences of all the genomes are concatenated, named `<genome>:<sequence id>`;
2. they are aligned against themselves with a single BLASTN;
3. the alignments of at least `--min_len` bp, not included in another alignment, are pooled by pair of genomes: the
   distance of two genomes is their number of errors (mismatches, plus one per extra fragment) per aligned bp, the worst
   of the two directions (`blast_distance_graph()` of `R/lib_blast.R`). The unaligned parts of the genomes are not
   counted.

# Installation

The pipeline depends on nextflow, that can be installed with:

```bash
curl -s https://get.nextflow.io | bash
```

# Usage

```bash
nextflow run pradosj/amr-tools/wf/evo/wg-blastn-distance --fasta='genomes/*.fasta'
```

Each FASTA file is a genome (a plasmid, a bacterial genome, ... with one or several sequences), named by its file name
(without the FASTA extension). A path suffixed by `!` is a FASTA file to split: each of its sequences is a genome,
named `<file name>:<sequence id>` (e.g. `--fasta='plasmids.fasta!'`). The options of BLASTN are given with
`--blastn_args` (default `-perc_identity 70 -evalue 1e-3 -subject_besthit -max_target_seqs 1000000`), the minimum
length of the alignments with `--min_len` (default 1000).

Outputs, in the output folder (`results/` by default):

- `wg-blastn-distance/report.html`: the HTML report, the minimum spanning tree of the distances (the genomes at distance
  0 merged into a single node);
- `wg-blastn-distance/distance_graph.rds`: the distance graph of the genomes (tbl_graph: an edge per pair of genomes
  aligning, with the number of mismatches, gap openings, fragments, aligned bp, errors and the distance (`weight`)).

# Test the pipeline

On 4 plasmids downloaded from ENA:

```bash
nextflow run pradosj/amr-tools/wf/evo/wg-blastn-distance/main.nf -profile standard,test
```
