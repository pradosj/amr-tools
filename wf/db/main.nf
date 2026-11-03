nextflow.enable.types = true

include { DATABASES; DB_ARCHIVE } from './subworkflows/db/main.nf'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'


// Build the databases of the pipelines in the output folder (checkm2/, org_db/, amrfinderplus/), and a tar archive of them
params {
	db: String? = null                  // databases folder (or archive) to complete: only the missing databases are built
	archive: String? = 'amr-tools-db'   // name of the archive (<archive>.tar.gz), null for no archive
}

workflow {
	main:
	// Validate parameters and print summary of supplied ones
	validateParameters()
	log.info(paramsSummaryLog(workflow))

	def dbs = DATABASES(params.db, null, null, null, true, true, true)
	def archive_ch = channel.empty()
	if (params.archive) {
		archive_ch = channel.empty().mix(DB_ARCHIVE(dbs.checkm2, dbs.org_db, dbs.amrfinderplus, params.archive))
	}

	publish:
	db = dbs.built
	archive = archive_ch
}

output {
	db: Channel<Record> {
		mode "copy"
		path { r -> r.dir >> "${r.name}" }
	}
	archive: Channel<Path> {
		mode "copy"
		path { x -> x >> "${x.name}" }
	}
}
