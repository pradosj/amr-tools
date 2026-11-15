# genome-rotate

`genome-rotate` a nextflow pipeline rotating the sequences of genomes (circular plasmids, chromosomes, ...) to start
at their best hit on a reference dataset (`GENOME_ROTATE`, see [`genome/rotate/blastn`](../../../registry/workflows/genome/rotate/blastn)):

1. a BLAST database is built from the reference dataset (e.g. the genes of the origins of replication, a reference
   plasmid);
2. each genome is aligned on it with BLASTN (the genome as query);
3. each sequence is rotated to start at its best hit (highest bit score), on the strand of the hit, and
   reverse-complemented when the hit is on the minus strand. The sequences without hit are left unchanged.

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
`--ref` is the FASTA file of the reference dataset (path or URL, optionally gzipped). The options of BLASTN are given
with `--blastn_args` (default `-perc_identity 70 -evalue 1e-3 -subject_besthit`).

Outputs, in the output folder (`results/` by default):

- `genome-rotate/<genome>.fasta`: the rotated genome;
- `genome-rotate/hits/<genome>.blastn`: the BLASTN hits of the genome (query) on the reference dataset (subject),
  `-outfmt "6 std qlen slen"`.

# Test the pipeline

Three plasmids rotated on a fourth one, downloaded from ENA:

```bash
nextflow run pradosj/amr-tools/wf/genome/rotate/main.nf -profile standard,test
```
