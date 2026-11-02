# assemblies-qc

`assemblies-qc` a nextflow pipeline for the quality control of bacterial assemblies: assembly statistics, alignment of the reads on the assembly, and completeness with CheckM2.

# Installation

The pipeline depends on nextflow, that can be installed with: 

```bash
curl -s https://get.nextflow.io | bash
```

# Usage

```bash
nextflow run pradosj/amr-tools/pipelines/assemblies-qc --fasta='assemblies/*.fasta' --ont='ont/*.fastq.gz' --illumina='illumina/*_R{1,2}.fastq.gz'
```

The sample id is taken from the file names: the reads of a sample are aligned on the assembly with the same sample id. Reads are optional.
Options: `--assemblies_qc.checkm2.skip`, `--assemblies_qc.checkm2.db` (CheckM2 database folder, downloaded when not given), `--assemblies_qc.checkm2.args`.

Outputs: `assemblies/<sample_id>__none.qc/checkm2_quality.tsv`, `reports/assemblies_multiqc.{txt,xlsx}` and `db/checkm2` when the database is downloaded.

# Test the pipeline

```bash
nextflow run pradosj/amr-tools/pipelines/assemblies-qc/main.nf -profile standard,test
```
