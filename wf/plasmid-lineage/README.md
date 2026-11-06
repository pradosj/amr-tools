# plasmid-lineage

`plasmid-lineage` a nextflow pipeline to study the lineage of a set of related plasmids:

1. the plasmids are aligned against each other (BLASTN) and rotated to start at the longest region (≥1 kb) shared by at least 80% of them, on its strand (`ROTATE_GENOMES`, see [`rotate-genomes`](../../workflows/rotate-genomes)); the next steps run on the rotated plasmids (`--no_rotation`: on the given plasmids);
2. the plasmids are annotated with `AMR_ANNOT` (ResFinder, PlasmidFinder, AMRFinderPlus, MOB-typer, ... see [`amr-annot`](../amr-annot); without orgfinder, and so without cgemlst) and prokka;
3. core genes: the prokka genes are clustered by reciprocal best BLASTN hits, and the clusters found in most plasmids are typed with pyMLST (allele distances, minimum spanning tree);
4. core genome: the regions of the plasmids shared by most plasmids;
5. an HTML report gathers the alignments, the annotations, the pyMLST typing and the core genome.

# Installation

The pipeline depends on nextflow, that can be installed with:

```bash
curl -s https://get.nextflow.io | bash
```

# Usage

```bash
nextflow run pradosj/amr-tools/wf/plasmid-lineage --fasta='plasmids/*.fasta' --name=my_plasmids --lineage.core_genome_count=3
```

Each FASTA file is a plasmid, named by its file name (without the FASTA extension), and all the plasmids form one set, named by `--name` (the sequence ids must be unique across the files). A gene or a region is core when it is shared by at least `--lineage.core_genome_count` plasmids (required, e.g. `--lineage.core_genome_count=3`). The options of the tools are given with `--lineage.blastn_args` (plasmids), `--lineage.genes_blastn_args` (genes, e.g. `-task blastn` for divergent genes), `--lineage.prokka_args` and `--lineage.pymlst_args`, and the AMR annotation with `--amr.<tool>.skip` / `--amr.<tool>.args`. The AMRFinderPlus database is read from the databases folder `--db.basedir` (default `results/db`), and downloaded when missing, see [Databases](../../README.md#databases).

Outputs: `<name>/` (`report.html`, the rotated plasmids `rotated/<plasmid>.fasta` (not with `--no_rotation`), the prokka annotation of each plasmid `prokka/<plasmid>`, the core genes tables, with `core_genes_index.tsv` (the plasmid of each index naming the plasmids in the core genes and pyMLST results), the core genome table, the pyMLST database and results), `amr/` and `reports/` (results and multi-sample report of `AMR_ANNOT`), `db/` (the databases built by the run).

# Test the pipeline

```bash
nextflow run pradosj/amr-tools/wf/plasmid-lineage/main.nf -profile standard,test
```
