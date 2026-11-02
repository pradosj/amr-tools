
# orgfinder

`orgfinder` a nextflow pipeline to detect organism in given assembled genome.

The pipeline uses FastANI to compare the assemblies to a reference database.

# Installation

The pipeline depends on nextflow, that can be installed with: 

```bash
curl -s https://get.nextflow.io | bash
```

# Usage

```bash
nextflow run pradosj/amr-tools/wf/orgfinder --query=query_assembly.fasta
```

# Test the pipeline

```bash
nextflow run pradosj/amr-tools/wf/orgfinder/main.nf -profile standard,test
```

