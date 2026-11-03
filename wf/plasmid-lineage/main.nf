nextflow.enable.types = true

include { PLASMID_LINEAGE; PlasmidLineageParams } from './workflows/plasmid-lineage/main.nf'
include { AmrParams } from './workflows/plasmid-lineage/workflows/amr-annot/main.nf'
include { DATABASES; DbParams } from './subworkflows/db/main.nf'
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

	def plasmids = channel.of(record(meta: [sample_id: params.name], fasta: files(params.fasta).toSorted()))
	def dbs = DATABASES(params.db, false, !params.amr?.orgfinder?.skip, !params.amr?.amrfinderplus?.skip)
	def lineage = PLASMID_LINEAGE(params.lineage, params.amr, dbs.org_db, dbs.amrfinderplus, plasmids)

	publish:
	rotated = lineage.rotated
	prokka = lineage.prokka
	core_genes = lineage.core_genes
	pymlst = lineage.pymlst
	core_genome = lineage.core_genome
	fastani = lineage.fastani
	report = lineage.report
	amr_multireport = lineage.amr_multireport
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
