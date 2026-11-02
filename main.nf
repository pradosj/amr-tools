nextflow.enable.types = true

include { SAMTOOLS_FASTQ as CONVERT_LONG_BAM_TO_FASTQ } from './modules/samtools/fastq'
include { READS_FILTER; ReadsFilterParams             } from './pipelines/reads-filter/main.nf'
include { READS_QC as FILTERED_READS_QC               } from './pipelines/reads-qc/main.nf'
include { ASSEMBLE; AssemblerParams                   } from './pipelines/assemble/main.nf'
include { ASSEMBLIES_QC; AssembliesQcParams           } from './pipelines/assemblies-qc/main.nf'
include { AMR_ANNOT; AmrParams                        } from './wf/amr_annot/main.nf'
include { IDENTITY as COPY_ASSEMBLY_FASTA             } from './modules/identity/main.nf'
include { IDENTITY as COPY_SHORT_READS                } from './modules/identity/main.nf'
include { IDENTITY as COPY_LONG_READS                 } from './modules/identity/main.nf'

// ------------------------------------------------------------------
// Main entry point when running the pipeline from command line
// ------------------------------------------------------------------
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

// Input reads given as glob patterns, the sample_id is taken from the file names
record ReadsParams {
	ont: String?
	illumina: String?
}

// Input assemblies given as glob patterns, the sample_id is taken from the file names
record AssembliesParams {
	fasta: String?
	info: String?
}

params {
	verbose: Boolean = false
	csv:     Path? = null
	reads:   ReadsParams = record(ont: null, illumina: null)
	assemblies: AssembliesParams = record(fasta: null, info: null)
	no_reads_filter: Boolean = false   // work with the given reads
	no_reads_qc: Boolean = false
	reads_filter: ReadsFilterParams = record(long_limit: 1_000_000_000, short_limit: 1_000_000_000, chopper_args: '--quality 15 --minlength 200')
	assembler: AssemblerParams = record(name: 'long_flye_medaka', args: record(flye: null, medaka: '--bacteria', spades: null, unicycler: null, hybracter: null))
	assemblies_qc: AssembliesQcParams = record(checkm2: record(skip: false, db: null, args: ''))
	amr: AmrParams = record(
		orgfinder:     record(skip: false, args: ''),
		speciator:     record(skip: true, args: ''),
		amrfinderplus: record(skip: false, args: ''),
		mobtyper:      record(skip: false, args: ''),
		plasmidfinder: record(skip: false, skip_longread: false, skip_shortread: false, args: ''),
		resfinder:     record(skip: false, skip_longread: false, skip_shortread: false, args: ''),
		cgemlst:       record(skip: false, args: ''),
		MLST:          record(skip: false, args: ''),    // args '' to autodetect the species
		prokka:        record(skip: true, args: '--kingdom Bacteria')
	)
}


// Input files of a sample
record Sample {
	sample_id: String
	reads_long: Path?
	reads_short_1: Path?
	reads_short_2: Path?
	assembly_fasta: Path?
}

def samplePart(id: String, reads_long: Path?, reads_short_1: Path?, reads_short_2: Path?, assembly_fasta: Path?) -> Sample {
	return record(sample_id: id, reads_long: reads_long, reads_short_1: reads_short_1, reads_short_2: reads_short_2, assembly_fasta: assembly_fasta)
}

// The file given for a field of a sample: each input file must be given only once
def uniqueFile(id: String, field: String, files: List<Path>) -> Path? {
	if (files.size() > 1) {
		def msg = "Sample ${id}: several files given for ${field}: ${files}"
		log.error(msg)    // errors raised in operators are only reported in .nextflow.log
		error(msg)
	}
	return files.isEmpty() ? null : files.first()
}

// Merge the parts of a sample found in the different input glob patterns
def mergeSample(id: String, parts: Bag<Sample>) -> Sample {
	return record(
		sample_id: id,
		reads_long: uniqueFile(id, 'reads_long', parts.collect({ p -> p.reads_long }).findAll({ f -> f }).toList()),
		reads_short_1: uniqueFile(id, 'reads_short_1', parts.collect({ p -> p.reads_short_1 }).findAll({ f -> f }).toList()),
		reads_short_2: uniqueFile(id, 'reads_short_2', parts.collect({ p -> p.reads_short_2 }).findAll({ f -> f }).toList()),
		assembly_fasta: uniqueFile(id, 'assembly_fasta', parts.collect({ p -> p.assembly_fasta }).findAll({ f -> f }).toList())
	)
}

// Sample from a row of the samplesheet (see assets/schema/samples_sheets_schema.json)
def sampleFromRow(row: Map) -> Sample {
	return samplePart(
		row.sample_id as String,
		row.reads_long ? file(row.reads_long as String) : null,
		row.reads_short_1 ? file(row.reads_short_1 as String) : null,
		row.reads_short_2 ? file(row.reads_short_2 as String) : null,
		row.assembly_fasta ? file(row.assembly_fasta as String) : null
	)
}


// Samples given with --csv, or found from the --reads.ont, --reads.illumina and --assemblies.fasta glob patterns
workflow SAMPLES {
	take:
	csv: Path?
	reads: ReadsParams
	assemblies: AssembliesParams

	main:
	def samples_ch = channel.empty()
	if (csv) {
		samples_ch = channel.fromList(samplesheetToList(csv, "${projectDir}/assets/schema/samples_sheets_schema.json"))
			.map({ row -> sampleFromRow(row[0] as Map) })
	} else {
		def illumina_re = '_(R?)([12])(_[0-9][0-9][0-9])?\\.(fq|fastq)\\.gz$'
		def parts_ch = channel.empty()
		if (reads.ont) {
			parts_ch = parts_ch.mix(channel.fromPath(reads.ont)
				.map({ f -> samplePart(f.name.replaceAll('\\.(fastq\\.gz|fq\\.gz|bam|cram)$', ''), f, null, null, null) }))
		}
		if (reads.illumina) {
			parts_ch = parts_ch.mix(channel.fromPath(reads.illumina)
				.map({ f ->
					def id = f.name.replaceAll(illumina_re, '')
					(f.name =~ illumina_re)[0][2] == '1' ? samplePart(id, null, f, null, null) : samplePart(id, null, null, f, null)
				}))
		}
		if (assemblies.fasta) {
			parts_ch = parts_ch.mix(channel.fromPath(assemblies.fasta)
				.map({ f -> samplePart(f.name.replaceAll('\\.(fasta|fa|fna)$', ''), null, null, null, f) }))
		}
		samples_ch = parts_ch
			.map({ s -> tuple(s.sample_id, s) })
			.groupBy()
			.map({ id, parts -> mergeSample(id, parts) })
	}

	emit:
	samples: Channel<Sample> = samples_ch
}


workflow {
	main:
	// Validate parameters and print summary of supplied ones
	validateParameters()
	log.info(paramsSummaryLog(workflow))

	def samples = SAMPLES(params.csv, params.reads, params.assemblies)

	// Copy the input reads in the work directory, so they can be published when not filtered
	def sr_rec = COPY_SHORT_READS(
		samples
			.filter({ s -> s.reads_short_1 })
			.map({ s -> record(meta: [sample_id: s.sample_id], file1: s.reads_short_1, file2: s.reads_short_2) })
	).map({ r -> record(meta: r.meta, reads: r.file1, read2: r.file2) })
	def lr_copy = COPY_LONG_READS(
		samples
			.filter({ s -> s.reads_long })
			.map({ s -> record(meta: [sample_id: s.sample_id], file1: s.reads_long, file2: null) })
	)

	// Convert long reads given in BAM/CRAM format into FASTQ format
	def bam_re = '.*\\.(bam|cram)$'
	def lr_rec = lr_copy
		.filter({ r -> !(r.file1.name ==~ bam_re) })
		.map({ r -> record(meta: r.meta, reads: r.file1) })
		.mix(CONVERT_LONG_BAM_TO_FASTQ(
			lr_copy
				.filter({ r -> r.file1.name ==~ bam_re })
				.map({ r -> record(meta: r.meta, bam: r.file1) })
		))

	// Filter reads
	if (!params.no_reads_filter) {
		def rd_flt_ch = READS_FILTER(params.reads_filter, sr_rec, lr_rec)
		sr_rec = rd_flt_ch.short_filtered
		lr_rec = rd_flt_ch.long_filtered
	}
	def reads_multiqc_ch = channel.empty()
	if (!params.no_reads_qc) {
		reads_multiqc_ch = FILTERED_READS_QC(sr_rec, lr_rec).multiqc
	}

	// De novo assembly of the samples given without assembly
	def todo_ch = samples
		.filter({ s -> !s.assembly_fasta })
		.map({ s -> record(meta: [sample_id: s.sample_id]) })
	def assemblies = ASSEMBLE(
		params.assembler,
		sr_rec.join(todo_ch, by: 'meta'),
		lr_rec.join(todo_ch, by: 'meta')
	)
	def given_asm_ch = COPY_ASSEMBLY_FASTA(
		samples
			.filter({ s -> s.assembly_fasta })
			.map({ s -> record(meta: [sample_id: s.sample_id], file1: s.assembly_fasta, file2: null) })
	)

	// Assemblies are identified by a meta map holding both sample_id and assembler_name
	def asm_rec = given_asm_ch
		.map({ r -> record(meta: r.meta + [assembler_name: 'none'], fasta: r.file1) })
		.mix(assemblies.map({ r -> record(meta: r.meta, fasta: r.fasta) }))
	def asm_dir_ch = assemblies.map({ r -> record(meta: r.meta, dir: r.dir) })

	// Reads of the sample of each assembly (a sample has a single assembly)
	def asm_keys = asm_rec.map({ r -> record(meta: r.meta.subMap(['sample_id']), asm_meta: r.meta) })
	def asm_short_reads = asm_keys.join(sr_rec, by: 'meta').map({ r -> record(meta: r.asm_meta, reads: r.reads, read2: r.read2) })
	def asm_long_reads = asm_keys.join(lr_rec, by: 'meta').map({ r -> record(meta: r.asm_meta, reads: r.reads) })

	// Compute assemblies QC stats
	def asm_qc_ch = ASSEMBLIES_QC(params.assemblies_qc, asm_rec, asm_short_reads, asm_long_reads)

	// Finally run AMR annotations
	def amr_report_ch = AMR_ANNOT(params.amr, asm_rec, sr_rec, lr_rec)

	publish:
	checkm2_db = asm_qc_ch.checkm2_db

	reads_long_filtered = lr_rec
	reads_short_filtered = sr_rec
	filtered_reads_multiqc = reads_multiqc_ch

	assemblies_fasta        = asm_rec
	assemblies_dir          = asm_dir_ch
	assemblies_checkm2      = asm_qc_ch.checkm2
	assemblies_multiqc      = asm_qc_ch.multiqc

	amr_multireport = amr_report_ch
}


output {
	checkm2_db: Channel<Path> {
		path { x -> x >> "db/checkm2" }
	}

	reads_long_filtered: Channel<Record> {
		path { r -> r.reads >> "reads/filtered_long/${r.meta.sample_id}.fastq.gz" }
	}
	reads_short_filtered: Channel<Record> {
		path { r ->
			r.reads >> "reads/filtered_short/${r.meta.sample_id}_R1.fastq.gz"
			r.read2 >> "reads/filtered_short/${r.meta.sample_id}_R2.fastq.gz"
		}
	}

	filtered_reads_multiqc: Channel<Record> {
		path { r -> r.html >> "reports/multiqc.html" }
	}

	assemblies_fasta: Channel<Record> {
		mode "copy" //because link is problematic
		path { r -> r.fasta >> "assemblies/${r.meta.sample_id}__${r.meta.assembler_name}.fasta" }
	}
	assemblies_dir: Channel<Record> {
		mode "copy" //because link is problematic
		enabled params.verbose
		path { r -> r.dir >> "assemblies/${r.meta.sample_id}__${r.meta.assembler_name}.output" }
	}
	assemblies_checkm2: Channel<Record> {
		path { r -> r.report >> "assemblies/${r.meta.sample_id}__${r.meta.assembler_name}.qc/checkm2_quality.tsv" }
	}
	assemblies_multiqc: Channel<Record> {
		path { r ->
			r.txt >> "reports/assemblies_multiqc.txt"
			r.xlsx >> "reports/assemblies_multiqc.xlsx"
		}
	}

	amr_multireport: Channel<Record> {
		mode "copy"    // the folder holds symlinks to the task outputs
		path { r ->
			r.folder.resolve('amr') >> "amr"    // results of the AMR tools, as gathered for the multireport
			r.html >> "reports/report.html"
			r.xlsx >> "reports/multitable.xlsx"
			r.assemblies_tsv >> "reports/multitable_assemblies.tsv"
			r.contigs_tsv >> "reports/multitable_contigs.tsv"
			r.resistances_tsv >> "reports/multitable_resistances.tsv"
			r.plasmidfinder_hits_tsv >> "reports/multitable_plasmidfinder_hits.tsv"
		}
	}
}
