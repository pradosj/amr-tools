
# gtf2gff

`gtf2gff` a helper tool to generate a GFF compatible with `bcftools` from the GTF of a bacterial genome obtained from NCBI.

# Installation

The pipeline depends on nextflow, that can be installed with: 

```bash
curl -s https://get.nextflow.io | bash
```

# Usage

```bash
nextflow run pradosj/amr-tools/pipelines/gtf2gff/main.nf --gtf 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/014/365/GCF_000014365.2_ASM1436v2/GCF_000014365.2_ASM1436v2_genomic.gtf.gz'
```

# Test the pipeline

```bash
nextflow run pradosj/amr-tools/pipelines/asm-compare/main.nf -profile standard,test
```

