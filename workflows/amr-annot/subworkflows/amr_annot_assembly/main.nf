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


// Options of a tool (params.amr.<tool>). Fields are nullable: modules fall back to their own defaults.
record ToolParams {
	skip: Boolean?
	args: String?
}

// Options of a tool with a database (params.amr.orgfinder, params.amr.amrfinderplus): db is the database folder
// (read by the pipelines, see DATABASES), args are given to the tool (fastANI for orgfinder)
record DbToolParams {
	skip: Boolean?
	db: Path?
	args: String?
}

// Options of AMR_ANNOT_ASSEMBLY (the assembly part of params.amr)
record AmrAssemblyParams {
	orgfinder: DbToolParams?
	speciator: ToolParams?
	amrfinderplus: DbToolParams?
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

	emit:
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
