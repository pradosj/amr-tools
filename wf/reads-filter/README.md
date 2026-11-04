# reads-filter

`reads-filter` a nextflow pipeline to filter sequencing reads before assembly.

ONT reads are filtered by quality and length with chopper. ONT and Illumina reads are then subsampled to a maximum number of bases.

# Installation

The pipeline depends on nextflow, that can be installed with: 

```bash
curl -s https://get.nextflow.io | bash
```

# Usage

```bash
nextflow run pradosj/amr-tools/wf/reads-filter --reads.ont='ont/*.fastq.gz' --reads.illumina='illumina/*_R{1,2}.fastq.gz'
```

ONT reads can also be given as BAM or CRAM files. The sample id is taken from the file names.
Options: `--reads_filter.long_limit` and `--reads_filter.short_limit` (maximum number of bases kept, negative to keep all reads), `--reads_filter.chopper_args` (default `--quality 15 --minlength 200`).

Outputs: `reads/filtered_long/<sample_id>.fastq.gz` and `reads/filtered_short/<sample_id>_R{1,2}.fastq.gz`.

# Test the pipeline

```bash
nextflow run pradosj/amr-tools/wf/reads-filter/main.nf -profile standard,test
```

The `test` profile uses ONT reads of *Klebsiella pneumoniae* downloaded from ENA ([ERR14747186](https://www.ebi.ac.uk/ena/browser/view/ERR14747186), 332 MB).
