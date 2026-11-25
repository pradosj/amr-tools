# assemblies-qc

`assemblies-qc` a nextflow pipeline for the quality control of bacterial assemblies: assembly statistics, alignment of the reads on the assembly, and completeness with CheckM2.

# Installation

The pipeline depends on nextflow, that can be installed with: 

```bash
curl -s https://get.nextflow.io | bash
```

# Usage

```bash
nextflow run pradosj/amr-tools/wf/genome/qc --assemblies.fasta='assemblies/*.fasta' --reads.ont='ont/*.fastq.gz' --reads.illumina='illumina/*_R{1,2}.fastq.gz'
```

The sample id is taken from the file names: the reads of a sample are aligned on the assembly with the same sample id. Reads are optional.
The topology of the contigs (circular or linear) given by the assemblers can be added with `--assemblies.json` (the `.asm.json` files published in `assemblies/` by the assembly pipelines, e.g. `--assemblies.fasta='results/assemblies/*.fasta' --assemblies.json='results/assemblies/*.asm.json'`): a JSON file goes with the FASTA file of the same sample id (the file name without `.asm.json` or `.json`); it is shown in the dashboard (topology of each contig, number of circular contigs).
Options: `--genome_qc.checkm2.skip`, `--genome_qc.checkm2.args`, and `--db.basedir` / `--db.checkm2` (databases folder, see [Databases](../../README.md#databases)). The CheckM2 database is downloaded when missing from the databases folder (`results/db` by default).

Outputs: `assemblies/<sample_id>__none.qc/checkm2_quality.tsv`, `assemblies/<sample_id>__none.qc/qc_report.html` (dashboard of the assembly: statistics, reads alignment, depth and topology of each contig), `reports/assemblies_multiqc.{txt,xlsx}` and `db/checkm2` when the database is downloaded (found by the next runs in the default databases folder `results/db`).

# Test the pipeline

```bash
nextflow run pradosj/amr-tools/wf/genome/qc/main.nf -profile standard,test
```
