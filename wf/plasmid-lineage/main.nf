nextflow.enable.types = true

include { PLASMID_LINEAGE; PlasmidLineageParams; PlasmidSet; CONCAT_FASTA; FIND_ALIGN_SEQ } from './workflows/plasmid-lineage/main.nf'
include { ROTATE_GENOMES } from './workflows/rotate-genomes/main.nf'
include { BLASTN as BLASTN_SELF } from './modules/blast/blastn'
include { AmrParams } from './workflows/plasmid-lineage/workflows/amr-annot/main.nf'
include { DATABASES; DbParams; Db } from './subworkflows/db/main.nf'
include { INPUT_FASTA } from './subworkflows/input/main.nf'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'


// FASTA files of the plasmids to compare (all the sequences form one set), their rotation, PLASMID_LINEAGE and AMR_ANNOT options
params {
	db: DbParams = record(basedir: 'results/db', amrfinderplus: 'amrfinderplus')    // databases folder (built when missing), or tar archive of it (path or URL)
	fasta: String
	name: String = 'plasmids'
	no_rotation: Boolean = false    // do not rotate the plasmids before the analysis
	lineage: PlasmidLineageParams = record(
		core_genome_count: null,    // required: --lineage.core_genome_count=<n>
		blastn_args: '-perc_identity 70 -evalue 1e-3 -subject_besthit',
		genes_blastn_args: '-subject_besthit',
		prokka_args: '--kingdom Bacteria',
		pymlst_args: '--coverage 0.85 --identity 0.85'
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

	// one set of all the plasmids (gzipped files are decompressed)
	def plasmids: Channel<PlasmidSet> = INPUT_FASTA(params.fasta)
		.map({ r -> r.fasta as Path })
		.collect()
		.flatMap({ fs -> fs.isEmpty() ? [] : [record(meta: record(sample_id: params.name), fasta: fs.toList().toSorted({ f -> f.name }))] })
	// Rotate the plasmids to start on the longest region shared by most of them, on the strand of this region
	// (unless --no_rotation)
	def lineage_input_ch = plasmids
	def rotated_ch = channel.empty()
	if (!params.no_rotation) {
		def blastn_args = params.lineage?.blastn_args ?: '-perc_identity 70 -evalue 1e-3 -subject_besthit'
		def fasta_ch = CONCAT_FASTA(plasmids)
		def self_ch = BLASTN_SELF(fasta_ch.map({ r -> record(query: record(meta: r.meta, path: r.fasta), subject: null, args: blastn_args) }))
			.map({ r -> record(meta: r.meta_query, tsv: r.tsv) })
		def align_seq_ch = FIND_ALIGN_SEQ(
			fasta_ch.join(self_ch, by: 'meta').map({ r -> record(meta: r.meta, fasta: r.fasta, blast: r.tsv) }),
			file("${projectDir}/workflows/plasmid-lineage/assets")
		)
		rotated_ch = ROTATE_GENOMES(
			fasta_ch,
			align_seq_ch.map({ r -> record(meta: r.meta, fasta: r.fasta, blast: null) }),
			blastn_args
		).rotated
		lineage_input_ch = rotated_ch.map({ r -> record(meta: r.meta, fasta: [r.fasta]) })
	}

	def dbs = DATABASES(params.db, (params.amr?.amrfinderplus?.skip ? [] : [Db.AMRFINDERPLUS]).toSet())    // no orgfinder in PLASMID_LINEAGE
	def lineage = PLASMID_LINEAGE(params.lineage, params.amr, dbs.amrfinderplus, lineage_input_ch)

	publish:
	rotated = rotated_ch
	prokka = lineage.prokka
	core_genes = lineage.core_genes
	pymlst = lineage.pymlst
	core_genome = lineage.core_genome
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
		path { r -> r.fasta >> "${r.meta.sample_id}/rotated.fasta" }
	}
	prokka: Channel<Record> {
		path { r -> r.dir >> "${r.meta.set}/prokka/${r.meta.sample_id}" }    // prokka of each plasmid
	}
	core_genes: Value<Record> {
		path { r ->
			r.genes_tsv >> "${r.meta.sample_id}/core_genes.tsv"
			r.clusters_tsv >> "${r.meta.sample_id}/core_clusters.tsv"
			r.index_tsv >> "${r.meta.sample_id}/core_genes_index.tsv"    // the plasmid of each index (pyMLST strain)
		}
	}
	pymlst: Value<Record> {
		path { r ->
			r.db >> "${r.meta.sample_id}/pymlst/pymlst.db"
			r.dir >> "${r.meta.sample_id}/pymlst/results"
		}
	}
	core_genome: Channel<Record> {
		path { r -> r.tsv >> "${r.meta.sample_id}/core_genome.tsv" }
	}
	report: Channel<Record> {
		path { r -> r.html >> "${r.meta.sample_id}/report.html" }
	}
	amr_assemblies: Channel<Record> {
		mode "copy"    // the folder holds symlinks to the task outputs
		path { r ->
			r.fasta >> "amr/assemblies/${r.meta.sample_id}__${r.meta.assembler_name}.fasta"
			r.fai >> "amr/assemblies/${r.meta.sample_id}__${r.meta.assembler_name}.fasta.fai"
			r.dir >> "amr/assemblies/${r.meta.sample_id}__${r.meta.assembler_name}.amr"    // results of the AMR tools on the assembly
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
