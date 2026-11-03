# plasmid-lineage

`plasmid-lineage` a nextflow pipeline to study the lineage of a set of related plasmids:

1. the plasmids are aligned against each other (BLASTN) and rotated to start at the longest region shared by at least 80% of them;
2. the rotated plasmids are annotated with `AMR_ANNOT` (ResFinder, PlasmidFinder, AMRFinderPlus, MOB-typer, ... see [`amr-annot`](../amr-annot)) and prokka;
3. core genes: the prokka proteins are clustered by reciprocal best BLASTP hits, and the clusters found in most plasmids are typed with pyMLST (allele distances, minimum spanning tree);
4. core genome: the regions of the rotated plasmids shared by most plasmids are compared with FastANI;
5. an HTML report gathers the alignments, the annotations, the pyMLST typing and the core genome ANI.

# Installation

The pipeline depends on nextflow, that can be installed with:

```bash
curl -s https://get.nextflow.io | bash
```

# Usage

```bash
nextflow run pradosj/amr-tools/wf/plasmid-lineage --fasta='plasmids/*.fasta' --name=my_plasmids
```

All the sequences of the FASTA files form one set of plasmids (sequence ids must be unique). A gene or a region is core when it is shared by at least `--lineage.core_genome_count` plasmids (default: 80% of the plasmids). The options of the tools are given with `--lineage.blastn_args`, `--lineage.blastp_args`, `--lineage.prokka_args`, `--lineage.pymlst_args` and `--lineage.fastani_args`, and the AMR annotation with `--amr.<tool>.skip` / `--amr.<tool>.args`.

Outputs: `<name>/` (`report.html`, rotated plasmids, prokka annotation, core genes and core genome tables, pyMLST database and results, FastANI of the core genomes), `amr/` and `reports/` (results and multi-sample report of `AMR_ANNOT`).

# Test the pipeline

```bash
nextflow run pradosj/amr-tools/wf/plasmid-lineage/main.nf -profile standard,test
```
