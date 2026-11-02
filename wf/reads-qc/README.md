# reads-qc

`reads-qc` a nextflow pipeline to check the quality of sequencing reads.

ONT reads are analysed with NanoPlot, Illumina reads with fastp, and a MultiQC report gathers all the samples.

# Installation

The pipeline depends on nextflow, that can be installed with: 

```bash
curl -s https://get.nextflow.io | bash
```

# Usage

```bash
nextflow run pradosj/amr-tools/wf/reads-qc --ont='ont/*.fastq.gz' --illumina='illumina/*_R{1,2}.fastq.gz'
```

The sample id is taken from the file names. Outputs: `<sample_id>/nanoplot/`, `<sample_id>/fastp.{json,html}` and `multiqc.html`.
