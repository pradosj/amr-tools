nextflow.enable.types = true

include { AMR_ANNOT_ASSEMBLY } from './subworkflows/amr_annot_assembly'
include { AMR_ANNOT_READS    } from './subworkflows/amr_annot_reads'
include { SAMTOOLS_FAIDX     } from './modules/samtools/faidx'
include { MULTIREPORT        } from './subworkflows/multireport'


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


workflow AMR_ANNOT {
	take:
	opts: AmrParams
	assemblies: Channel<Genome>       // meta: [sample_id, assembler_name]
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
