nextflow.enable.types = true

include { PLASMID_LINEAGE; PlasmidLineageParams } from './workflows/plasmid-lineage/main.nf'
include { AmrParams } from './workflows/plasmid-lineage/workflows/amr-annot/main.nf'
include { DATABASES; DbParams; Db } from './subworkflows/db/main.nf'
include { INPUT_FASTA } from './subworkflows/input/main.nf'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'


// FASTA files of the plasmids to compare (glob, all the sequences form one set), PLASMID_LINEAGE and AMR_ANNOT options
params {
	db: DbParams = record(basedir: 'results/db', checkm2: 'checkm2', org_db: 'org_db', amrfinderplus: 'amrfinderplus')    // databases folder (built when missing), or tar archive of it (path or URL)
	fasta: String
	name: String = 'plasmids'
	lineage: PlasmidLineageParams = record(
		core_genome_count: null,    // null: 80% of the plasmids
		blastn_args: '-perc_identity 70 -evalue 1e-3 -subject_besthit',
		blastp_args: '-subject_besthit',
		prokka_args: '--kingdom Bacteria',
		pymlst_args: '--coverage 0.85 --identity 0.85',
		fastani_args: '--fragLen 500'
	)
	amr: AmrParams = record(
		orgfinder:     record(skip: true, args: ''),
		speciator:     record(skip: true, args: ''),
		amrfinderplus: record(skip: false, args: ''),
		mobtyper:      record(skip: false, args: ''),
		plasmidfinder: record(skip: false, skip_longread: true, skip_shortread: true, args: ''),
		resfinder:     record(skip: false, skip_longread: true, skip_shortread: true, args: ''),
		cgemlst:       record(skip: true, args: ''),
		MLST:          record(skip: true, args: ''),
		prokka:        record(skip: true, args: '')    // prokka is run by PLASMID_LINEAGE
	)
}

workflow {
	main:
	// Validate parameters and print summary of supplied ones
	validateParameters()
	log.info(paramsSummaryLog(workflow))

	// one set of all the plasmids (gzipped files are decompressed)
	def plasmids = INPUT_FASTA(params.fasta)
		.map({ r -> r.fasta as Path })
		.collect()
		.flatMap({ fs -> fs.isEmpty() ? [] : [record(meta: [sample_id: params.name], fasta: fs.toList().toSorted({ f -> f.name }))] })
	def dbs = DATABASES(params.db, ((params.amr?.orgfinder?.skip ? [] : [Db.ORG_DB]) + (params.amr?.amrfinderplus?.skip ? [] : [Db.AMRFINDERPLUS])).toSet())
	def lineage = PLASMID_LINEAGE(params.lineage, params.amr, dbs.org_db, dbs.amrfinderplus, plasmids)

	publish:
	rotated = lineage.rotated
	prokka = lineage.prokka
	core_genes = lineage.core_genes
	pymlst = lineage.pymlst
	core_genome = lineage.core_genome
	fastani = lineage.fastani
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
		path { r -> r.dir >> "${r.meta.sample_id}/prokka" }
	}
	core_genes: Channel<Record> {
		path { r ->
			r.genes_tsv >> "${r.meta.sample_id}/core_genes.tsv"
			r.clusters_tsv >> "${r.meta.sample_id}/core_clusters.tsv"
		}
	}
	pymlst: Channel<Record> {
		path { r ->
			r.db >> "${r.meta.sample_id}/pymlst/pymlst.db"
			r.dir >> "${r.meta.sample_id}/pymlst/results"
		}
	}
	core_genome: Channel<Record> {
		path { r -> r.tsv >> "${r.meta.sample_id}/core_genome.tsv" }
	}
	fastani: Channel<Record> {
		path { r -> r.tsv >> "${r.meta.sample_id}/core_genome_fastani.tsv" }
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
