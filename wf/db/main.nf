nextflow.enable.types = true

include { DATABASES; DB_FOLDER; DbParams; Db } from './subworkflows/db/main.nf'
include { TAR                  } from './modules/tar/tar.nf'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'


// Build the databases of the pipelines in the db/ output folder (checkm2/, org_db/, amrfinderplus/), and a tar archive of them
params {
	db: DbParams = record(basedir: 'results/db', checkm2: 'checkm2', org_db: 'org_db', amrfinderplus: 'amrfinderplus')    // databases folder (or archive) to complete: only the missing databases are built
	archive: String? = 'amr-tools-db'   // name of the archive (<archive>.tar.gz), null for no archive
}

workflow {
	main:
	validateParameters()
	log.info(paramsSummaryLog(workflow))

	def dbs = DATABASES(params.db, [Db.CHECKM2, Db.ORG_DB, Db.AMRFINDERPLUS].toSet())
	def archive_ch = channel.empty()
	if (params.archive) {
		def folder = DB_FOLDER(dbs.checkm2, dbs.org_db, dbs.amrfinderplus)
		archive_ch = channel.empty().mix(TAR(folder.map({ d -> record(dir: d, archive_name: "${params.archive}.tar.gz") })))
	}

	publish:
	db = dbs.built
	archive = archive_ch
}

output {
	db: Channel<Record> {
		mode "copy"
		path { r -> r.dir >> "db/${r.name}" }    // the default databases folder of the pipelines (--db.basedir=results/db)
	}
	archive: Channel<Path> {
		mode "copy"
		path { x -> x >> "${x.name}" }
	}
}
