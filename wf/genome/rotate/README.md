# genome-rotate

`genome-rotate` a nextflow pipeline rotating the sequences of genomes (circular plasmids, chromosomes, ...) to start
at their alignment on a reference dataset (`GENOME_ROTATE`, see [`genome/rotate`](../../../registry/workflows/genome/rotate)):

1. each genome is aligned with minimap2 (`-x asm5`) on the reference dataset (e.g. the genes of the origins of
   replication, a reference plasmid);
2. each sequence is rotated to start at its longest alignment on the reference sequence it covers most, on the strand
   of the alignment, and reverse-complemented when it is on the minus strand. The sequences without alignment are
   left unchanged.

# Installation

The pipeline depends on nextflow, that can be installed with:

```bash
curl -s https://get.nextflow.io | bash
```

# Usage

```bash
nextflow run pradosj/amr-tools/wf/genome/rotate --fasta='genomes/*.fasta' --ref=reference.fasta
```

Each FASTA file is a genome, named by its file name (without the FASTA extension). A path suffixed by `!` is a FASTA
file to split: each of its sequences is a genome, named `<file name>:<sequence id>` (e.g. `--fasta='plasmids.fasta!'`).
`--ref` is the FASTA file of the reference dataset (path or URL, optionally gzipped).


Outputs, in the output folder (`results/` by default):

- `genome-rotate/<genome>.fasta`: the rotated genome;
- `genome-rotate/rotation/<genome>.tsv`: the rotation of each sequence of the genome (columns `chr`, `pos`, `strand`;
  `pos` 1 and strand `+` for the sequences left unchanged).

# Test the pipeline

Three plasmids rotated on a fourth one, downloaded from ENA:

```bash
nextflow run pradosj/amr-tools/wf/genome/rotate/main.nf -profile standard,test
```
