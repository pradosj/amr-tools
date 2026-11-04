nextflow.enable.types = true

include { RESFINDER             } from './modules/cgetools/resfinder'
include { PLASMIDFINDER         } from './modules/cgetools/plasmidfinder'
include { ORGFINDER; OrgfinderDb  } from './workflows/orgfinder/main.nf'
include { SPECIATOR             } from './modules/speciator'
include { AMRFINDERPLUS_RUN     } from './modules/amrfinderplus/run'
include { PROKKA                } from './modules/tseemann/prokka'
include { MLST                  } from './modules/tseemann/mlst'
include { MLST as CGEMLST       } from './modules/cgetools/mlst'
include { MOBTYPER              } from './modules/mobsuite/mobtyper.nf'
include { SAMTOOLS_FAIDX        } from './modules/samtools/faidx'


// Options of a tool (params.amr.<tool>). Fields are nullable: modules fall back to their own defaults.
record ToolParams {
	skip: Boolean?
	args: String?
}

// Options of AMR_ANNOT_ASSEMBLY (the assembly part of params.amr)
record AmrAssemblyParams {
	orgfinder: ToolParams?
	speciator: ToolParams?
	amrfinderplus: ToolParams?
	mobtyper: ToolParams?
	plasmidfinder: ToolParams?
	resfinder: ToolParams?
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


// Gather the AMR results of an assembly into a folder (published as amr/assemblies/<sample_id>__<assembler_name>.amr), read into an RDS file
process AMR_COMBINE_ASSEMBLY {
	container 'registry.gitlab.unige.ch/amr-genomics/rscript:main'
	label 'rscript'
	memory 8.GB
	cpus 2
	time 30.m
	input:
	record(meta: Map, fasta: Path, org_name: String?, orgfinder: Path?, amrfinderplus: Path?, resfinder: Path?, plasmidfinder: Path?, cgemlst: Path?, mobtyper: Path?, MLST: Path?, prokka: Path?)
	libs: List<Path>    // lib_io.R, lib_typing.R
	stage:
	stageAs(libs, 'lib/*')
	stageAs(fasta, 'assembly.fasta')
	stageAs(orgfinder, 'amr/orgfinder')
	stageAs(amrfinderplus, 'amr/amrfinderplus')
	stageAs(resfinder, 'amr/resfinder')
	stageAs(plasmidfinder, 'amr/plasmidfinder')
	stageAs(cgemlst, 'amr/mlst')
	stageAs(mobtyper, 'amr/mobtyper.tsv')
	stageAs(MLST, 'amr/MLST.txt')
	stageAs(prokka, 'amr/prokka')
	output:
	record(meta: meta, dir: file('amr'), rds: file('amr.rds'))
	script:
	"""
	#!/usr/bin/env Rscript
	source("lib/lib_io.R")
	source("lib/lib_typing.R")
	fs::dir_create("amr")
	${org_name != null ? "writeLines('${org_name}','amr/org_name.txt')" : ''}
	read_assemblies_amr("amr","assembly.fasta",'${meta.sample_id}','${meta.assembler_name}') |> saveRDS("amr.rds")
	"""
	stub:
	"""
	mkdir -p amr
	touch amr.rds
	"""
}


workflow AMR_ANNOT_ASSEMBLY {
	take:
	opts: AmrAssemblyParams
	org_db: Value<Path?>              // orgfinder database folder (null when orgfinder is skipped)
	amrfinderplus_db: Value<Path?>    // AMRFinderPlus database folder (null when AMRFinderPlus is skipped)
	assemblies: Channel<Genome>

	main:
	def org_args = new groovy.json.JsonSlurper().parseText(file("${moduleDir}/assets/default_org_args.json").text) as Map

	// CGE - RESFINDER
	def resfinder_ch = channel.empty()
	if (!opts.resfinder?.skip) {
		resfinder_ch = RESFINDER(assemblies.map({ r -> r + record(args: opts.resfinder?.args) }))
	}

	// Plasmid typing
	def plasmidfinder_ch = channel.empty()
	if (!opts.plasmidfinder?.skip) {
		plasmidfinder_ch = PLASMIDFINDER(assemblies.map({ r -> r + record(args: opts.plasmidfinder?.args) }))
	}

	// NCBI AMRfinder+
	def amrfinderplus_ch = channel.empty()
	if (!opts.amrfinderplus?.skip) {
		amrfinderplus_ch = AMRFINDERPLUS_RUN(
			assemblies.combine(db: amrfinderplus_db).map({ r -> r + record(args: opts.amrfinderplus?.args) })
		)
	}

	// MOBsuite - MOBtyper
	def mobtyper_ch = channel.empty()
	if (!opts.mobtyper?.skip) {
		mobtyper_ch = MOBTYPER(assemblies.map({ r -> r + record(args: opts.mobtyper?.args) }))
	}

	// PROKKA annotations (disabled unless skip is explicitly false)
	def prokka_ch = channel.empty()
	if (opts.prokka?.skip == false) {
		prokka_ch = PROKKA(assemblies.map({ r -> r + record(args: opts.prokka?.args) }))
	}

	// Speciator (disabled unless skip is explicitly false)
	def speciator_ch = channel.empty()
	if (opts.speciator?.skip == false) {
		speciator_ch = SPECIATOR(assemblies)
	}

	def MLST_ch = channel.empty()
	if (!opts.MLST?.skip) {
		MLST_ch = MLST(assemblies.map({ r -> r + record(args: opts.MLST?.args) }))
	}

	// Run orgfinder to auto detect organism
	def orgfinder_ch = channel.empty()
	if (!opts.orgfinder?.skip) {
		def db: Value<OrgfinderDb> = org_db.map({ d -> record(dir: d, fasta: files("${d}/fna/*.fna").toSorted(), tsv: (d as Path).resolve('db.tsv')) })
		orgfinder_ch = ORGFINDER(db, assemblies, opts.orgfinder?.args).orgfinder
	}

	// MLST typing on detected org_name, when there is a scheme for it
	def cgemlst_ch = channel.empty()
	if (!opts.cgemlst?.skip) {
		cgemlst_ch = CGEMLST(
			assemblies
				.join(orgfinder_ch.map({ r -> record(meta: r.meta, org_name: r.org_name) }), by: 'meta')
				.map({ r -> record(
					meta: r.meta,
					fasta: r.fasta,
					species_args: (org_args[r.org_name] as Map)?.cgemlst as String,
					args: opts.cgemlst?.args
				) })
				.filter({ r -> r.species_args })
		)
	}

	// Gather the results of each assembly into a folder and an RDS file
	def fai_ch = SAMTOOLS_FAIDX(assemblies)
	def combine_ch = AMR_COMBINE_ASSEMBLY(
		assemblies
			.join(orgfinder_ch.map({ r -> record(meta: r.meta, org_name: r.org_name, orgfinder: r.dir) }), by: 'meta', remainder: true)
			.join(amrfinderplus_ch.map({ r -> record(meta: r.meta, amrfinderplus: r.dir) }), by: 'meta', remainder: true)
			.join(resfinder_ch.map({ r -> record(meta: r.meta, resfinder: r.dir) }), by: 'meta', remainder: true)
			.join(plasmidfinder_ch.map({ r -> record(meta: r.meta, plasmidfinder: r.dir) }), by: 'meta', remainder: true)
			.join(cgemlst_ch.map({ r -> record(meta: r.meta, cgemlst: r.dir) }), by: 'meta', remainder: true)
			.join(mobtyper_ch.map({ r -> record(meta: r.meta, mobtyper: r.tsv) }), by: 'meta', remainder: true)
			.join(MLST_ch.map({ r -> record(meta: r.meta, MLST: r.txt) }), by: 'meta', remainder: true)
			.join(prokka_ch.map({ r -> record(meta: r.meta, prokka: r.dir) }), by: 'meta', remainder: true),
		[file("${moduleDir}/assets/lib_io.R"), file("${moduleDir}/assets/lib_typing.R")]
	)
	def combined_ch = assemblies
		.map({ r -> record(meta: r.meta, fasta: r.fasta) })
		.join(fai_ch.map({ r -> record(meta: r.meta, fai: r.fai) }), by: 'meta')
		.join(combine_ch, by: 'meta')

	emit:
	combined: Channel<Record>      = combined_ch        // records with fields meta, fasta, fai, dir (results folder), rds (read_assemblies_amr() table)
	orgfinder: Channel<Record>     = orgfinder_ch       // records with fields meta, dir, org_name, org_acc, ani
	amrfinderplus: Channel<Record> = amrfinderplus_ch   // records with fields meta, dir
	resfinder: Channel<Record>     = resfinder_ch       // records with fields meta, dir
	mobtyper: Channel<Record>      = mobtyper_ch        // records with fields meta, tsv
	plasmidfinder: Channel<Record> = plasmidfinder_ch   // records with fields meta, dir
	prokka: Channel<Record>        = prokka_ch          // records with fields meta, dir
	speciator: Channel<Record>     = speciator_ch       // records with fields meta, dir
	cgemlst: Channel<Record>       = cgemlst_ch         // records with fields meta, dir
	MLST: Channel<Record>          = MLST_ch            // records with fields meta, txt
}
