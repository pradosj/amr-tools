
# asm-compare

`asm-compare` a simple nextflow pipeline to compare two assemblies and generate a VCF.

The pipeline steps include:

  1) Circular rotation of the query sequences to match with best reference sequence hit (this step can be skip with option --no_query_rotation)
  
  2) Map query assemblies to the reference genome with `minimap2 -x asm5` and produce an indexed `bam` file

  3) Generate a VCF listing the mutations with `bcftools mpileup`
  
  4) Generate a simple 3 columns table list of mutations from the VCF


# Installation

The pipeline depends on nextflow, that can be installed with: 

```bash
curl -s https://get.nextflow.io | bash
```

# Usage

```bash
nextflow run pradosj/amr-tools/pipelines/asm-compare/main.nf --ref=ref_assembly.fasta --query=query_assembly.fasta
```

# Test the pipeline

```bash
nextflow run pradosj/amr-tools/pipelines/asm-compare/main.nf -profile standard,test
```

