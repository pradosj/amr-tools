# assemble

`assemble` a nextflow pipeline for the de novo assembly of bacterial genomes from ONT and/or Illumina reads.

# Installation

The pipeline depends on nextflow, that can be installed with: 

```bash
curl -s https://get.nextflow.io | bash
```

# Usage

```bash
nextflow run pradosj/amr-tools/wf/assemble --reads.ont='ont/*.fastq.gz' --assembler.name=long_flye_medaka
nextflow run pradosj/amr-tools/wf/assemble --reads.ont='ont/*.fastq.gz' --reads.illumina='illumina/*_R{1,2}.fastq.gz' --assembler.name=hybrid_unicycler
```

ONT reads can also be given as BAM or CRAM files. The sample id is taken from the file names.
Assembly strategies (`--assembler.name`): `long_flye`, `long_flye_medaka` (default), `long_unicycler`, `long_hybracter`, `short_unicycler`, `short_spades`, `hybrid_unicycler`, `hybrid_hybracter`, `hybrid_flye_medaka_pilon`. Tool arguments: `--assembler.args.<tool>`.

Outputs: `assemblies/<sample_id>__<assembler>.fasta`, and with `--verbose` the assembler output folder `assemblies/<sample_id>__<assembler>.output`.

# Test the pipeline

```bash
nextflow run pradosj/amr-tools/wf/assemble/main.nf -profile standard,test
```

The `test` profile uses ONT reads of *Klebsiella pneumoniae* downloaded from ENA ([ERR14747186](https://www.ebi.ac.uk/ena/browser/view/ERR14747186), 332 MB).
