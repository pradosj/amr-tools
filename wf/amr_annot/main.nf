nextflow.enable.types = true

include { RESFINDER                                } from './modules/cgetools/resfinder'
include { RESFINDER     as RESFINDER_LONGREAD      } from './modules/cgetools/resfinder'
include { RESFINDER     as RESFINDER_SHORTREAD     } from './modules/cgetools/resfinder'
include { PLASMIDFINDER                            } from './modules/cgetools/plasmidfinder'
include { PLASMIDFINDER as PLASMIDFINDER_LONGREAD  } from './modules/cgetools/plasmidfinder'
include { PLASMIDFINDER as PLASMIDFINDER_SHORTREAD } from './modules/cgetools/plasmidfinder'
include { ORGFINDER_DETECT                         } from './modules/orgfinder/detect'
include { SPECIATOR                                } from './modules/speciator'
include { AMRFINDERPLUS_UPDATE                     } from './modules/amrfinderplus/update'
include { AMRFINDERPLUS_RUN                        } from './modules/amrfinderplus/run'
include { PROKKA                                   } from './modules/tseemann/prokka'
include { MLST                                     } from './modules/tseemann/mlst'
include { MLST as CGEMLST                          } from './modules/cgetools/mlst'
include { MOBTYPER                                 } from './modules/mobsuite/mobtyper.nf'
include { SAMTOOLS_FAIDX                           } from './modules/samtools/faidx'
include { MULTIREPORT; NamedEntry                  } from './subworkflows/multireport'


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
record Assembly {
	meta: Map
	fasta: Path
}

// Short reads sample: read2 is null for single-end reads
record ShortReads {
	meta: Map
	reads: Path
	read2: Path?
}

// Long reads sample: one FASTQ per sample
record LongReads {
	meta: Map
	reads: Path
}


// Entries of the AMR results folder read by the multireport (see lib_typing.R)
def assemblyEntry(meta: Map, suffix: String, path: Path) -> NamedEntry {
	return record(name: "amr/assemblies/${meta.sample_id}__${meta.assembler_name}${suffix}", path: path, text: null)
}

def amrEntry(meta: Map, path: Path) -> NamedEntry {
	return assemblyEntry(meta, ".amr/${path.name}", path)
}

def readsEntry(meta: Map, name: String, path: Path) -> NamedEntry {
	return record(name: "amr/reads/${meta.sample_id}.reads_amr/${name}", path: path, text: null)
}

// A result of an AMR tool, published under the given tool name
record AmrResult {
	meta: Map
	tool: String
	path: Path
}

def amrResult(meta: Map, tool: String, path: Path) -> AmrResult {
	return record(meta: meta, tool: tool, path: path)
}


workflow AMR_ANNOT_ASSEMBLY {
	take:
	opts: AmrParams
	assemblies: Channel<Assembly>

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
		def amrfinderplus_db = AMRFINDERPLUS_UPDATE()
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
		orgfinder_ch = ORGFINDER_DETECT(assemblies.map({ r -> r + record(args: opts.orgfinder?.args) }))
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


workflow AMR_ANNOT_READS {
	take:
	opts: AmrParams
	short_reads: Channel<ShortReads>
	long_reads: Channel<LongReads>

	main:
	def plasmidfinder_long_ch = channel.empty()
	if (!opts.plasmidfinder?.skip_longread) {
		plasmidfinder_long_ch = PLASMIDFINDER_LONGREAD(
			long_reads.map({ r -> record(meta: r.meta, long_reads: r.reads, args: opts.plasmidfinder?.args) })
		)
	}
	def resfinder_long_ch = channel.empty()
	if (!opts.resfinder?.skip_longread) {
		resfinder_long_ch = RESFINDER_LONGREAD(
			long_reads.map({ r -> record(meta: r.meta, long_reads: r.reads, args: opts.resfinder?.args) })
		)
	}
	def plasmidfinder_short_ch = channel.empty()
	if (!opts.plasmidfinder?.skip_shortread) {
		plasmidfinder_short_ch = PLASMIDFINDER_SHORTREAD(short_reads.map({ r -> r + record(args: opts.plasmidfinder?.args) }))
	}
	def resfinder_short_ch = channel.empty()
	if (!opts.resfinder?.skip_shortread) {
		resfinder_short_ch = RESFINDER_SHORTREAD(short_reads.map({ r -> r + record(args: opts.resfinder?.args) }))
	}

	emit:
	resfinder_long: Channel<Record>      = resfinder_long_ch        // records with fields meta, dir
	resfinder_short: Channel<Record>     = resfinder_short_ch       // records with fields meta, dir
	plasmidfinder_long: Channel<Record>  = plasmidfinder_long_ch    // records with fields meta, dir
	plasmidfinder_short: Channel<Record> = plasmidfinder_short_ch   // records with fields meta, dir
}


workflow AMR_ANNOT {
	take:
	opts: AmrParams
	assemblies: Channel<Assembly>       // meta: [sample_id, assembler_name]
	short_reads: Channel<ShortReads>    // meta: [sample_id]
	long_reads: Channel<LongReads>      // meta: [sample_id]

	main:
	def fai_ch = SAMTOOLS_FAIDX(assemblies)
	def asm = AMR_ANNOT_ASSEMBLY(opts, assemblies)
	def rds = AMR_ANNOT_READS(opts, short_reads, long_reads)

	// Gather all the results in a folder for the multireport
	def entries = assemblies.map({ r -> assemblyEntry(r.meta, '.fasta', r.fasta) })
		.mix(fai_ch.map({ r -> assemblyEntry(r.meta, '.fasta.fai', r.fai) }))
		.mix(asm.orgfinder.map({ r -> record(name: "amr/assemblies/${r.meta.sample_id}__${r.meta.assembler_name}.amr/org_name.txt", path: null, text: r.org_name as String) }))
		.mix(asm.orgfinder.map({ r -> amrEntry(r.meta, r.dir) }))
		.mix(asm.amrfinderplus.map({ r -> amrEntry(r.meta, r.dir) }))
		.mix(asm.resfinder.map({ r -> amrEntry(r.meta, r.dir) }))
		.mix(asm.plasmidfinder.map({ r -> amrEntry(r.meta, r.dir) }))
		.mix(asm.cgemlst.map({ r -> amrEntry(r.meta, r.dir) }))
		.mix(asm.mobtyper.map({ r -> amrEntry(r.meta, r.tsv) }))
		.mix(asm.MLST.map({ r -> amrEntry(r.meta, r.txt) }))
		.mix(asm.prokka.map({ r -> amrEntry(r.meta, r.dir) }))
		.mix(rds.resfinder_long.map({ r -> readsEntry(r.meta, 'resfinder_long', r.dir) }))
		.mix(rds.plasmidfinder_long.map({ r -> readsEntry(r.meta, 'plasmidfinder_long', r.dir) }))
		.mix(rds.resfinder_short.map({ r -> readsEntry(r.meta, 'resfinder_short', r.dir) }))
		.mix(rds.plasmidfinder_short.map({ r -> readsEntry(r.meta, 'plasmidfinder_short', r.dir) }))
	def report = MULTIREPORT(entries)

	// Results of the tools, by assembly and by reads sample
	def asm_results_ch = asm.orgfinder.map({ r -> amrResult(r.meta, 'orgfinder', r.dir) })
		.mix(asm.amrfinderplus.map({ r -> amrResult(r.meta, 'amrfinderplus', r.dir) }))
		.mix(asm.resfinder.map({ r -> amrResult(r.meta, 'resfinder', r.dir) }))
		.mix(asm.mobtyper.map({ r -> amrResult(r.meta, 'mobtyper', r.tsv) }))
		.mix(asm.plasmidfinder.map({ r -> amrResult(r.meta, 'plasmidfinder', r.dir) }))
		.mix(asm.cgemlst.map({ r -> amrResult(r.meta, 'cgemlst', r.dir) }))
		.mix(asm.MLST.map({ r -> amrResult(r.meta, 'MLST', r.txt) }))
		.mix(asm.prokka.map({ r -> amrResult(r.meta, 'prokka', r.dir) }))
	def reads_results_ch = rds.resfinder_long.map({ r -> amrResult(r.meta, 'resfinder_long', r.dir) })
		.mix(rds.resfinder_short.map({ r -> amrResult(r.meta, 'resfinder_short', r.dir) }))
		.mix(rds.plasmidfinder_long.map({ r -> amrResult(r.meta, 'plasmidfinder_long', r.dir) }))
		.mix(rds.plasmidfinder_short.map({ r -> amrResult(r.meta, 'plasmidfinder_short', r.dir) }))

	emit:
	assemblies_results: Channel<AmrResult> = asm_results_ch        // meta: [sample_id, assembler_name]
	reads_results: Channel<AmrResult>      = reads_results_ch      // meta: [sample_id]
	multireport: Channel<Record>           = report                // records with fields folder, html, xlsx, assemblies_tsv, contigs_tsv, resistances_tsv, plasmidfinder_hits_tsv
}
