nextflow.enable.types = true

include { PLASMID_LINEAGE; PlasmidLineageParams; Genome } from './workflows/plasmid-lineage/main.nf'
include { GENOME_ROTATE } from './workflows/rotate/main.nf'
include { ROTATION_ORIGINS } from './subworkflows/rotation_origins/main.nf'
include { IDENTITY as COPY_FASTA } from './modules/identity/main.nf'
include { AmrParams } from './workflows/plasmid-lineage/workflows/amr-annot/main.nf'
include { DATABASES; DbParams; Db } from './subworkflows/db/main.nf'
include { INPUT_FASTA; asIdPath } from './subworkflows/input/main.nf'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'


// FASTA files of the plasmids to compare (one plasmid per file, all forming the set --name), their rotation, PLASMID_LINEAGE
// and AMR_ANNOT options
params {
	db: DbParams = record(basedir: 'results/db', amrfinderplus: 'amrfinderplus')    // databases folder (built when missing), or tar archive of it (path or URL)
	fasta: String
	name: String = 'plasmids'
	no_rotation: Boolean = false    // do not rotate the plasmids before the analysis
	lineage: PlasmidLineageParams = record(
		blastn_args: '-perc_identity 70 -evalue 1e-3 -subject_besthit'
	)
	amr: AmrParams = record(
		amrfinderplus: record(skip: false, args: ''),
		mobtyper:      record(skip: false, args: ''),
		plasmidfinder: record(skip: false, args: ''),
		resfinder:     record(skip: false, args: ''),
		prokka:        record(skip: true, args: '')
	)
}

workflow {
	main:
	// Validate parameters and print summary of supplied ones
	validateParameters()
	log.info(paramsSummaryLog(workflow))

	// the plasmids: one per FASTA file, named by its sample id (gzipped files are decompressed), copied in the work
	// directory (to be published when they are not rotated)
	def plasmids: Channel<Genome> = COPY_FASTA(INPUT_FASTA(params.fasta).map({ r -> record(meta: record(id: r.meta.id), file1: r.path, file2: null) }))
		.map({ r -> record(meta: r.meta, fasta: r.file1) })
	// Rotate each plasmid to start on the longest region shared by most of them, on the strand of this region (unless
	// --no_rotation)
	def lineage_input_ch = plasmids
	def rotated_ch = channel.empty()
	if (!params.no_rotation) {
		def origins_ch = ROTATION_ORIGINS(
			record(id: params.name),
			plasmids.map({ r -> asIdPath(record(meta: r.meta, path: r.fasta)) }),
			params.lineage?.blastn_args
		)
		// the same origins for each plasmid
		rotated_ch = GENOME_ROTATE(
			plasmids.combine(origins_ch).map({ p, a -> record(fasta: record(meta: p.meta, path: p.fasta), origins: a) })
		).fasta
		lineage_input_ch = rotated_ch.map({ r -> record(meta: r.meta, fasta: r.path) })
	}

	def dbs = DATABASES(params.db, (params.amr?.amrfinderplus?.skip ? [] : [Db.AMRFINDERPLUS]).toSet())    // no orgfinder in PLASMID_LINEAGE
	def lineage = PLASMID_LINEAGE(params.lineage, params.amr, dbs.amrfinderplus, params.name, lineage_input_ch)

	publish:
	rotated = rotated_ch
	report = lineage.report
	amr_assemblies = lineage.amr_assemblies
	amr_multiqc = lineage.amr_multiqc
	db = dbs.built
}

output {
	db: Channel<Record> {
		path { r -> r.dir >> "db/${r.name}" }    // databases built by this run: use --db=<output dir>/db in the next runs
	}
	rotated: Channel<Record> {
		path { r -> r.path >> "${params.name}/rotated/${r.meta.id}.fasta" }    // each rotated plasmid
	}
	report: Channel<Record> {
		path { r -> r.html >> "${r.meta.id}/report.html" }
	}
	amr_assemblies: Channel<Record> {
		mode "copy"    // the folder holds symlinks to the task outputs
		path { r ->
			r.fasta >> "amr/assemblies/${r.meta.id}__${r.meta.assembler_name}.fasta"
			r.fai >> "amr/assemblies/${r.meta.id}__${r.meta.assembler_name}.fasta.fai"
			r.dir >> "amr/assemblies/${r.meta.id}__${r.meta.assembler_name}.amr"    // results of the AMR tools on the assembly
		}
	}
	amr_multiqc: Channel<Record> {
		path { r ->
			r.html >> "reports/report.html"
			r.rds >> "reports/amr_multiqc.rds"    // AMR results of all the assemblies and samples, as read by lib_typing.R
			r.xlsx >> "reports/multitable.xlsx"
			r.assemblies_tsv >> "reports/multitable_assemblies.tsv"
			r.contigs_tsv >> "reports/multitable_contigs.tsv"
			r.resistances_tsv >> "reports/multitable_resistances.tsv"
			r.plasmidfinder_hits_tsv >> "reports/multitable_plasmidfinder_hits.tsv"
		}
	}
}
