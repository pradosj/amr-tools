# genome-qc

`genome-qc` a nextflow pipeline for the quality control of bacterial assemblies: assembly statistics, alignment of the reads on the assembly, and completeness with CheckM2.

# Installation

The pipeline depends on nextflow, that can be installed with: 

```bash
curl -s https://get.nextflow.io | bash
```

# Usage

```bash
nextflow run pradosj/amr-tools/wf/genome/qc --genomes.fasta='genomes/*.fasta' --reads.ont='ont/*.fastq.gz' --reads.illumina='illumina/*_R{1,2}.fastq.gz'
```

The sample id is taken from the file names: the reads of a sample are aligned on the assembly with the same sample id. Reads are optional.
The topology of the contigs (circular or linear) given by the assemblers can be added with `--genomes.json` (the `.info.json` files published in `genomes/` by the assembly pipelines, e.g. `--genomes.fasta='results/genomes/*.fasta' --genomes.json='results/genomes/*.info.json'`): a JSON file goes with the FASTA file of the same sample id (the file name without `.info.json` or `.json`); it is shown in the dashboard (topology of each contig, number of circular contigs).
Options: `--genome_qc.checkm2.skip`, `--genome_qc.checkm2.args`, and `--db.basedir` / `--db.checkm2` (databases folder, see [Databases](../../README.md#databases)). The CheckM2 database is downloaded when missing from the databases folder (`results/db` by default).

Outputs: `genomes/<sample_id>__none.qc/checkm2_quality.tsv`, `genomes/<sample_id>__none.qc/qc_report.html` (dashboard of the assembly: statistics, reads alignment, depth and topology of each contig), `reports/genomes_multiqc.{txt,xlsx}` and `db/checkm2` when the database is downloaded (found by the next runs in the default databases folder `results/db`).

# Test the pipeline

```bash
nextflow run pradosj/amr-tools/wf/genome/qc/main.nf -profile standard,test
```
