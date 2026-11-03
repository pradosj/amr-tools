# amr-annot

`amr-annot` a nextflow pipeline for the antimicrobial resistance and plasmid annotation of bacterial assemblies (orgfinder, AMRFinderPlus, ResFinder, PlasmidFinder, MOB-typer, MLST, ...) and of their reads (ResFinder, PlasmidFinder), with a multi-sample report.

# Installation

The pipeline depends on nextflow, that can be installed with: 

```bash
curl -s https://get.nextflow.io | bash
```

# Usage

```bash
nextflow run pradosj/amr-tools/wf/amr-annot --assemblies.fasta='assemblies/*.fasta' --reads.ont='ont/*.fastq.gz' --reads.illumina='illumina/*_R{1,2}.fastq.gz'
```

The sample id is taken from the file names. Reads are optional. Each tool can be enabled or disabled with `--amr.<tool>.skip` and given arguments with `--amr.<tool>.args`.

orgfinder compares the assemblies to a database of NCBI reference genomes, and AMRFinderPlus uses its own database. A database that is not given is built (orgfinder) or downloaded (AMRFinderPlus) and saved to `db/`; give them back with `--db=results/db` (or `--amr.orgfinder.db=results/db/org_db`, `--amr.amrfinderplus.db=results/db/amrfinderplus`) to skip the build, see [Databases](../../README.md#databases).

Outputs: `amr/` (results of the tools, by assembly and by reads sample), `reports/` (HTML report and tables) and `db/` (the databases built by the run).

# Test the pipeline

```bash
nextflow run pradosj/amr-tools/wf/amr-annot/main.nf -profile standard,test
```
