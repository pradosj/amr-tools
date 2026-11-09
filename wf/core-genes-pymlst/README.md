# core-genes-pymlst

`core-genes-pymlst` a nextflow pipeline to find the core genes of a set of genomes and type them with pyMLST
(`CORE_GENES_PYMLST`, see [`core-genes-pymlst`](../../workflows/core-genes-pymlst)):

1. the genes of each genome are annotated with prokka;
2. the genes of each pair of genomes are aligned with BLASTN (each genome against itself included, for the paralogs),
   and their reciprocal hits are kept (`BLAST_RECIPROCAL_HITS`, see
   [`blast_reciprocal_hits`](../../subworkflows/blast_reciprocal_hits));
3. the genes are clustered by their reciprocal hits, with one gene per genome in each cluster: the clusters found in
   at least `--n_core` genomes are the core genes;
4. the core genes are typed with pyMLST (wgMLST): allele distances between the genomes, MLST table, alignment;
5. an HTML report shows the distances, the alleles and the minimum spanning tree of the genomes.

# Installation

The pipeline depends on nextflow, that can be installed with:

```bash
curl -s https://get.nextflow.io | bash
```

# Usage

```bash
nextflow run pradosj/amr-tools/wf/core-genes-pymlst --fasta='genomes/*.fasta' --n_core=3
```

Each FASTA file is a genome (a plasmid, a bacterial genome, ... with one or several sequences), named by its file name
(without the FASTA extension). A path suffixed by `!` is a FASTA file to split: each of its sequences is a genome,
named `<file name>:<sequence id>` (e.g. `--fasta='plasmids.fasta!'`). A gene is core when it is found in at least `--n_core` genomes (required). The options
of the tools are given with `--prokka_args`, `--blastn_args` (e.g. `-task blastn` for divergent genes) and
`--pymlst_args`.

Outputs, in the output folder (`results/` by default):

- `core_genes_index.tsv`: the genome of each index: the genomes are named by their index (from 1, in the order of
  their meta) in the other outputs;
- `core_genes.tsv` (the genes of the clusters, `<index>:<prokka id>`, with their genome, prokka product, length,
  cluster, and whether the cluster is core) and `core_clusters.tsv` (the size and number of genomes of each cluster);
- `pymlst/`: the pyMLST database (`pymlst.db`), its template (`template.fasta`, the core genes of the genome with the
  most core genes) and results (`results/`: `distance.tsv`, `mlst.tsv`, `msa.fasta`), the strains being the indexes
  of the genomes;
- `prokka/<genome>`: the prokka annotation of each genome;
- `report.html`: a heatmap of the pyMLST distances, an interactive heatmap (plotly) of the alleles of the core genes
  in the genomes, and the minimum spanning tree of the genomes for the pyMLST distances.

# Test the pipeline

On 4 plasmids downloaded from ENA (`--n_core=4`):

```bash
nextflow run pradosj/amr-tools/wf/core-genes-pymlst/main.nf -profile standard,test
```
