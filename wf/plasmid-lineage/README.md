# plasmid-lineage

`plasmid-lineage` a nextflow pipeline to study the lineage of a set of related plasmids:

1. the plasmids are aligned against each other (BLASTN) and rotated to start at the longest region (≥1 kb) shared by at least 80% of them, on its strand (`ROTATE_GENOMES`, see [`rotate-genomes`](../../workflows/rotate-genomes)); the next steps run on the rotated plasmids (`--no_rotation`: on the given plasmids);
2. the plasmids are annotated with `AMR_ANNOT` (ResFinder, PlasmidFinder, AMRFinderPlus, MOB-typer, ... see [`amr-annot`](../amr-annot); without orgfinder, and so without cgemlst) and prokka;
3. core genes: the prokka proteins are clustered by reciprocal best BLASTP hits, and the clusters found in most plasmids are typed with pyMLST (allele distances, minimum spanning tree);
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

All the sequences of the FASTA files form one set of plasmids (sequence ids must be unique). A gene or a region is core when it is shared by at least `--lineage.core_genome_count` plasmids (required, e.g. `--lineage.core_genome_count=3`). The options of the tools are given with `--lineage.blastn_args`, `--lineage.blastp_args`, `--lineage.prokka_args` and `--lineage.pymlst_args`, and the AMR annotation with `--amr.<tool>.skip` / `--amr.<tool>.args`. The AMRFinderPlus database is read from the databases folder `--db.basedir` (default `results/db`), and downloaded when missing, see [Databases](../../README.md#databases).

Outputs: `<name>/` (`report.html`, rotated plasmids `rotated.fasta` (not with `--no_rotation`), prokka annotation, core genes and core genome tables, pyMLST database and results), `amr/` and `reports/` (results and multi-sample report of `AMR_ANNOT`), `db/` (the databases built by the run).

# Test the pipeline

```bash
nextflow run pradosj/amr-tools/wf/plasmid-lineage/main.nf -profile standard,test
```
