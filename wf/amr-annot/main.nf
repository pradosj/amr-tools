nextflow.enable.types = true

include { AMR_ANNOT_ASSEMBLY } from './subworkflows/amr_annot_assembly'
include { AMR_ANNOT_READS    } from './subworkflows/amr_annot_reads'
include { SAMTOOLS_FAIDX     } from './modules/samtools/faidx'
include { MULTIREPORT        } from './subworkflows/multireport'
include { READS_INPUT; ASSEMBLIES_INPUT } from './conf/inputs.nf'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'


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
	def report = MULTIREPORT(
		assemblies.map({ r -> record(name: "amr/assemblies/${r.meta.sample_id}__${r.meta.assembler_name}.fasta", path: r.fasta) })
			.mix(fai_ch.map({ r -> record(name: "amr/assemblies/${r.meta.sample_id}__${r.meta.assembler_name}.fasta.fai", path: r.fai) }))
			.mix(asm.orgfinder.map({ r -> record(name: "amr/assemblies/${r.meta.sample_id}__${r.meta.assembler_name}.amr/org_name.txt", path: null, text: r.org_name) }))
			.mix(asm.orgfinder.map({ r -> record(name: "amr/assemblies/${r.meta.sample_id}__${r.meta.assembler_name}.amr/${r.dir.name}", path: r.dir) }))
			.mix(asm.amrfinderplus.map({ r -> record(name: "amr/assemblies/${r.meta.sample_id}__${r.meta.assembler_name}.amr/${r.dir.name}", path: r.dir) }))
			.mix(asm.resfinder.map({ r -> record(name: "amr/assemblies/${r.meta.sample_id}__${r.meta.assembler_name}.amr/${r.dir.name}", path: r.dir) }))
			.mix(asm.plasmidfinder.map({ r -> record(name: "amr/assemblies/${r.meta.sample_id}__${r.meta.assembler_name}.amr/${r.dir.name}", path: r.dir) }))
			.mix(asm.cgemlst.map({ r -> record(name: "amr/assemblies/${r.meta.sample_id}__${r.meta.assembler_name}.amr/${r.dir.name}", path: r.dir) }))
			.mix(asm.mobtyper.map({ r -> record(name: "amr/assemblies/${r.meta.sample_id}__${r.meta.assembler_name}.amr/${r.tsv.name}", path: r.tsv) }))
			.mix(asm.MLST.map({ r -> record(name: "amr/assemblies/${r.meta.sample_id}__${r.meta.assembler_name}.amr/${r.txt.name}", path: r.txt) }))
			.mix(asm.prokka.map({ r -> record(name: "amr/assemblies/${r.meta.sample_id}__${r.meta.assembler_name}.amr/${r.dir.name}", path: r.dir) }))
			.mix(rds.resfinder_long.map({ r -> record(name: "amr/reads/${r.meta.sample_id}.reads_amr/resfinder_long", path: r.dir) }))
			.mix(rds.plasmidfinder_long.map({ r -> record(name: "amr/reads/${r.meta.sample_id}.reads_amr/plasmidfinder_long", path: r.dir) }))
			.mix(rds.resfinder_short.map({ r -> record(name: "amr/reads/${r.meta.sample_id}.reads_amr/resfinder_short", path: r.dir) }))
			.mix(rds.plasmidfinder_short.map({ r -> record(name: "amr/reads/${r.meta.sample_id}.reads_amr/plasmidfinder_short", path: r.dir) }))
	)

	emit:
	multireport: Channel<Record> = report                // records with fields folder, html, xlsx, assemblies_tsv, contigs_tsv, resistances_tsv, plasmidfinder_hits_tsv
}


// Input assemblies and reads given as glob patterns (the sample_id is taken from the file names), and AMR_ANNOT options
params {
	fasta: String? = null
	ont: String? = null
	illumina: String? = null
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

workflow {
	main:
	// Validate parameters and print summary of supplied ones
	validateParameters()
	log.info(paramsSummaryLog(workflow))

	def assemblies = ASSEMBLIES_INPUT(params.fasta)
	def reads = READS_INPUT(params.ont, params.illumina)
	def report = AMR_ANNOT(params.amr, assemblies, reads.short_reads, reads.long_reads)

	publish:
	amr_multireport = report
}

output {
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
