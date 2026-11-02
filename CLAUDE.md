# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

`amr-tools` is a Nextflow pipeline for processing sequencing data of antimicrobial multi-resistant bacteria, with a focus on plasmids. It runs reads filtering → reads QC → de novo assembly → assembly QC → AMR annotation. Every step is also a standalone pipeline in `pipelines/`. Everything is written in **typed syntax** (`nextflow.enable.types = true` in every script); requires Nextflow `>=26.10.4` and the `nf-schema` plugin.

## Commands

```bash
# Run on test assemblies (assets/tests/*.fasta) with docker
nextflow run main.nf -profile standard,test

# Run on inputs (params are nested records: --reads.ont, --reads.illumina, --assemblies.fasta)
nextflow run main.nf --reads.ont='reads/*.fastq.gz' --assembler.name=long_flye_medaka
nextflow run main.nf --csv samples.csv            # samplesheet, validated by assets/schema/samples_sheets_schema.json

# Standalone pipelines (flat --ont / --illumina / --fasta globs, same option names as main.nf)
nextflow run pipelines/amr-annot/main.nf -profile standard,test
nextflow run pipelines/reads-qc/main.nf --ont='ont/*.fastq.gz' --illumina='ill/*_R{1,2}.fastq.gz'

# Dry run using process `stub:` blocks (use -without-docker when docker cannot mount the work dir)
nextflow run main.nf -profile test -stub -without-docker

# Lint — keep it warning-free
nextflow lint .
```

`nextflow lint` is not enough: several typed-syntax errors only show up at run time, so always also run a `-stub` test after a change. Values starting with `--` must be given with `=` on the CLI (`--mpileup_args='--no-BAQ'`). `-entry` is not supported with the strict parser.

`resume = true` is set in `nextflow.config`, so runs resume by default. Outputs go to `results/` (change with `-output-dir`). Profiles (shared by all pipelines): `standard`/`docker` (docker, linux/amd64), `arm64` (arm64 medaka image), `hpc` (slurm + apptainer, UNIGE cluster paths), `tiny` (low resource limits), plus a pipeline specific `test`.

## Architecture

- `main.nf`: the entry workflow. It declares all its params with their defaults in a typed `params {}` block (nested params are `record` types exported by the pipelines, e.g. `AssemblerParams`, `AmrParams`), validates them with `nf-schema`, builds a channel of `Sample` records with the `SAMPLES` workflow (from `--csv`, or from the `reads`/`assemblies` glob params), converts ONT BAM/CRAM to FASTQ, then chains the named workflows `READS_FILTER`, `READS_QC`, `ASSEMBLE`, `ASSEMBLIES_QC` and `AMR_ANNOT` included from `pipelines/*/main.nf`. `--no_reads_filter` / `--no_reads_qc` skip the first two. It also declares the `publish:`/`output {}` blocks, so all result paths of the full pipeline are defined there.
- `pipelines/<name>/`: one folder per step (`reads-filter`, `reads-qc`, `assemble`, `assemblies-qc`, `amr-annot`) and per extra tool (`orgfinder`, `asm-compare`, `gtf2gff`). Each `main.nf` holds the named workflow(s) and their record types, plus its own `params {}`, entry workflow and `output {}` used only when run directly (an included file's params and entry workflow are ignored). Each folder has its own `nextflow.config`, `nextflow_schema.json`, `README.md`, `assets/` (scripts referenced via `${moduleDir}`) and maybe nested `subworkflows/`. Pipeline composition with `include { workflow as X }` is not used: a composed pipeline only takes plain values, not channels.
- `modules/<tool>/...`: the single canonical copy of each process. **Pipelines and subworkflows include modules through relative symlinks** (e.g. `pipelines/amr-annot/subworkflows/amr_annot_assembly/modules/cgetools -> ../../../../../modules/cgetools`). To use a new module, add a symlink under the including folder's `modules/` dir. `..` in include paths is resolved textually through symlinks, so always include from a local `./modules/...` link.
- `conf/` (reached from every pipeline through a `conf -> ../../conf` symlink): `profiles.config` (shared profiles), `containers.config` (shared containers, e.g. `withLabel: 'rscript'`), and `inputs.nf` (`READS_INPUT` and `ASSEMBLIES_INPUT`: reading of the glob inputs of the standalone pipelines, Illumina R1/R2 pairing, ONT BAM/CRAM conversion; its modules are under `conf/modules/`).
- `nextflow.config` files only hold setup, manifest, the `includeConfig` of `conf/` and the `test` profile. Retry/resource escalation is encoded in the processes (`errorStrategy task.attempt < 2 ? 'retry' : 'ignore'`, `memory 8.GB + 8.GB * (task.attempt - 1)`), not in `withName` selectors. Light tasks use `label 'local'` (run locally by the `hpc` profile), R scripts use `label 'rscript'`.
- Module-local executables live in `modules/<tool>/resources/usr/bin/` (enabled by `nextflow.enable.moduleBinaries`).

## Typed syntax conventions

- A module takes ONE destructured record input holding data and options, e.g. `record(meta: Map, reads: Path, read2: Path?, args: String?)`, and emits a `record(...)` (see `modules/fastp/main.nf`). Modules never read `meta`, they pass it through.
- Tool arguments are explicit `args: String?` inputs used as `${args ?: ''}` (no `task.ext.args`). The caller always passes them (from params, or `null`): typed processes have no input defaults and no omitted inputs.
- Records: `ShortReads { meta; reads: Path; read2: Path? }` (`read2` null for single-end), `LongReads { meta; reads: Path }`, `Assembly { meta; fasta: Path }`. Reads meta is `[sample_id: ...]`; assembly meta is `[sample_id: ..., assembler_name: ...]` (`assembler_name: 'none'` for assemblies given as input) and output paths use `${sample_id}__${assembler_name}`. Add fields with `r + record(args: ...)`, value inputs with `ch.combine(db: value)`; join with `join(other, by: 'meta')`. Record types are structural: a record with extra fields can be passed where a narrower record type is expected.
- Every field of a params record must be nullable and the code must fall back to defaults on null: a partial CLI override (`--amr.prokka.skip=false`) replaces the whole record.
- Casts (`as Map`, `as Path`) are only needed to read a property or use an operator on an untyped field (`?`), not to pass it to a function, record or process.
- Known pitfalls: `files(glob)` is unordered (use `.toSorted()`); typed staging stages the same source file only once even under two names; a workflow emit named like a local variable fails at run time; a workflow with a single emit (or a pipeline with a single publish) returns the value itself; a line starting with `+` is a new statement; `mix` of records with different field sets may need a `field: null`.

## Remaining untyped code

Unused modules not yet converted: `modules/rscript`, `modules/ncbi/taxdump`, `modules/dorado`, `modules/gzip`, `modules/kraken2`, `modules/fastqc`, `modules/pling`. Convert a module to typed syntax before using it.
