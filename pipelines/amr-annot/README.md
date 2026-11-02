# amr-annot

`amr-annot` a nextflow pipeline for the antimicrobial resistance and plasmid annotation of bacterial assemblies (orgfinder, AMRFinderPlus, ResFinder, PlasmidFinder, MOB-typer, MLST, ...) and of their reads (ResFinder, PlasmidFinder), with a multi-sample report.

# Installation

The pipeline depends on nextflow, that can be installed with: 

```bash
curl -s https://get.nextflow.io | bash
```

# Usage

```bash
nextflow run pradosj/amr-tools/pipelines/amr-annot --fasta='assemblies/*.fasta' --ont='ont/*.fastq.gz' --illumina='illumina/*_R{1,2}.fastq.gz'
```

The sample id is taken from the file names. Reads are optional. Each tool can be enabled or disabled with `--amr.<tool>.skip` and given arguments with `--amr.<tool>.args`.

Outputs: `amr/` (results of the tools, by assembly and by reads sample) and `reports/` (HTML report and tables).

# Test the pipeline

```bash
nextflow run pradosj/amr-tools/pipelines/amr-annot/main.nf -profile standard,test
```
