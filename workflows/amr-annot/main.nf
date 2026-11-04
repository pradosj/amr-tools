nextflow.enable.types = true

include { AMR_ANNOT_ASSEMBLY } from './subworkflows/amr_annot_assembly'
include { AMR_ANNOT_READS    } from './subworkflows/amr_annot_reads'
include { RMD_RENDER         } from './modules/rmd/render'


// Options of a tool (params.amr.<tool>). Fields are nullable: modules fall back to their own defaults.
record ToolParams {
	skip: Boolean?
	args: String?
}

// Options of a tool that can also run on the reads
record ReadsToolParams {
	skip: Boolean?             // skip on the assemblies
	skip_longread: Boolean?    // skip on the long reads
	skip_shortread: Boolean?   // skip on the short reads
	args: String?
}

// Options of AMR_ANNOT (params.amr)
record AmrParams {
	orgfinder: ToolParams?
	speciator: ToolParams?
	amrfinderplus: ToolParams?
	mobtyper: ToolParams?
	plasmidfinder: ReadsToolParams?
	resfinder: ReadsToolParams?
	cgemlst: ToolParams?
	MLST: ToolParams?
	prokka: ToolParams?
}

// Assembly to annotate. meta holds both sample_id and assembler_name.
record Genome {
	meta: Map
	fasta: Path
	gff: Path?
	info: Path?
}

// Short reads sample: read2 is null for single-end reads
record ShortReads {
	meta: Map?
	read1: Path
	read2: Path?
}

// Long reads sample: one FASTQ per sample
record LongReads {
	meta: Map?
	reads: Path
}


// Aggregate the AMR results of all the assemblies and reads samples: summary tables, and an RDS file of the results (for the report)
process AMR_MULTIQC {
	container 'registry.gitlab.unige.ch/amr-genomics/rscript:main'
	label 'rscript'
	memory 8.GB
	cpus 2
	time 30.m
	input:
	record(assemblies: List<Path>, reads: List<Path>)    // RDS files of AMR_COMBINE_ASSEMBLY and AMR_COMBINE_READS
	libs: List<Path>    // lib_io.R, lib_typing.R
	stage:
	stageAs(libs, 'lib/*')
	stageAs(assemblies, 'assemblies/amr*.rds')
	stageAs(reads, 'reads/amr*.rds')
	output:
	record(
		rds: file('amr_multiqc.rds'),
		xlsx: file('multitable.xlsx'),
		assemblies_tsv: file('multitable_assemblies.tsv'),
		contigs_tsv: file('multitable_contigs.tsv'),
		resistances_tsv: file('multitable_resistances.tsv'),
		plasmidfinder_hits_tsv: file('multitable_plasmidfinder_hits.tsv')
	)
	script:
	"""
	#!/usr/bin/env Rscript
	source("lib/lib_io.R")
	source("lib/lib_typing.R")
	db <- db_aggregate(Sys.glob("assemblies/amr*.rds"),Sys.glob("reads/amr*.rds"))
	saveRDS(db,"amr_multiqc.rds")
	tbl <- amr_tables(db)
	openxlsx::write.xlsx(tbl,file="multitable.xlsx")
	readr::write_tsv(tbl\$assemblies,file="multitable_assemblies.tsv")
	readr::write_tsv(tbl\$contigs,file="multitable_contigs.tsv")
	readr::write_tsv(tbl\$resistances,file="multitable_resistances.tsv")
	readr::write_tsv(tbl\$plasmidfinder_hits,file="multitable_plasmidfinder_hits.tsv")
	"""
	stub:
	"""
	touch amr_multiqc.rds multitable.xlsx multitable_assemblies.tsv multitable_contigs.tsv multitable_resistances.tsv multitable_plasmidfinder_hits.tsv
	"""
}


workflow AMR_ANNOT {
	take:
	opts: AmrParams
	org_db: Value<Path?>              // orgfinder database folder (null when orgfinder is skipped)
	amrfinderplus_db: Value<Path?>    // AMRFinderPlus database folder (null when AMRFinderPlus is skipped)
	assemblies: Channel<Genome>       // meta: [sample_id, assembler_name]
	short_reads: Channel<ShortReads>    // meta: [sample_id]
	long_reads: Channel<LongReads>      // meta: [sample_id]

	main:
	def libs = [file("${moduleDir}/assets/lib_io.R"), file("${moduleDir}/assets/lib_typing.R")]
	def asm = AMR_ANNOT_ASSEMBLY(opts, org_db, amrfinderplus_db, assemblies)
	def rds = AMR_ANNOT_READS(opts, short_reads, long_reads)

	// Aggregate the results of all the assemblies and reads samples (no table without results)
	def multiqc_ch = AMR_MULTIQC(
		asm.combined.map({ r -> record(assembly: true, rds: r.rds) })
			.mix(rds.combined.map({ r -> record(assembly: false, rds: r.rds) }))
			.collect()
			.flatMap({ rs -> rs.isEmpty() ? [] : [record(
				assemblies: rs.findAll({ r -> r.assembly }).collect({ r -> r.rds as Path }).toSorted(),
				reads: rs.findAll({ r -> !r.assembly }).collect({ r -> r.rds as Path }).toSorted()
			)] }),
		libs
	)
	def html_ch = RMD_RENDER(multiqc_ch.map({ r -> record(
		template: file("${moduleDir}/assets/multireport.Rmd"),
		inputs: [r.rds] + libs,
		render_params: "db_rds='${r.rds.name}'"
	) }))

	emit:
	assemblies_amr: Channel<Record> = asm.combined                                       // records with fields meta, fasta, fai, dir, rds
	reads_amr: Channel<Record>      = rds.combined                                       // records with fields meta, dir, rds
	multiqc: Channel<Record>        = multiqc_ch.combine(html_ch).map({ t, h -> t + record(html: h) })    // records with fields rds, html, xlsx, assemblies_tsv, contigs_tsv, resistances_tsv, plasmidfinder_hits_tsv
}
