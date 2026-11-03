nextflow.enable.types = true

include { CHECKM2_DB           } from './modules/checkm2/db'
include { AMRFINDERPLUS_UPDATE } from './modules/amrfinderplus/update'
include { UNTAR                } from './modules/tar/untar.nf'
include { ORG_DB_BUILD         } from './subworkflows/org_db/main.nf'


// A database built by DATABASES, published under its name in the databases folder
record BuiltDb {
	name: String    // checkm2, org_db or amrfinderplus
	dir: Path
}


// Archive the databases folder, with checkm2/, org_db/ and amrfinderplus/ at its root
process DB_ARCHIVE {
	container 'docker.io/ubuntu:24.04'
	memory 4.GB
	cpus 1
	time 4.h

	input:
	checkm2: Path
	org_db: Path
	amrfinderplus: Path
	name: String

	stage:
	stageAs(checkm2, 'checkm2')
	stageAs(org_db, 'org_db')
	stageAs(amrfinderplus, 'amrfinderplus')

	output:
	file("${name}.tar.gz")

	script:
	"""
	tar -czhf '${name}.tar.gz' checkm2 org_db amrfinderplus
	"""

	stub:
	"""
	touch '${name}.tar.gz'
	"""
}


// Databases of the pipelines. Each database is, by priority: the given folder (checkm2_dir, org_db_dir,
// amrfinderplus_dir), its sub-folder in `db`, or built (downloaded). `db` is a databases folder (as built by
// wf/db, or the db/ folder of the results of a previous run) or a tar archive of it (path or URL), with the
// sub-folders checkm2/, org_db/ and amrfinderplus/; an archive must hold all the needed databases.
// A database that is not needed is null.
workflow DATABASES {
	take:
	db: String?
	checkm2_dir: Path?
	org_db_dir: Path?
	amrfinderplus_dir: Path?
	need_checkm2: Boolean
	need_org_db: Boolean
	need_amrfinderplus: Boolean

	main:
	def is_archive = db != null && db ==~ /.+\.(tar|tar\.gz|tgz|tar\.bz2|tar\.xz)$/
	def dir = (db != null && !is_archive) ? file(db) : null
	def root: Value<Path?> = channel.value(null)
	if (is_archive && (need_checkm2 || need_org_db || need_amrfinderplus)) {
		root = UNTAR(channel.value(record(archive: file(db), dir_name: 'db')))
	}
	def built_ch: Channel<BuiltDb> = channel.empty()

	// CheckM2
	def checkm2_ch: Value<Path?> = channel.value(null)
	if (need_checkm2) {
		if (checkm2_dir != null) {
			checkm2_ch = channel.value(checkm2_dir)
		} else if (is_archive) {
			checkm2_ch = root.map({ d -> (d as Path).resolve('checkm2') })
		} else if (dir != null && dir.resolve('checkm2').exists()) {
			checkm2_ch = channel.value(dir.resolve('checkm2'))
		} else {
			def checkm2_built = CHECKM2_DB()
			checkm2_ch = checkm2_built
			built_ch = built_ch.mix(checkm2_built.map({ d -> record(name: 'checkm2', dir: d) }))
		}
	}

	// orgfinder
	def org_db_ch: Value<Path?> = channel.value(null)
	if (need_org_db) {
		if (org_db_dir != null) {
			org_db_ch = channel.value(org_db_dir)
		} else if (is_archive) {
			org_db_ch = root.map({ d -> (d as Path).resolve('org_db') })
		} else if (dir != null && dir.resolve('org_db').exists()) {
			org_db_ch = channel.value(dir.resolve('org_db'))
		} else {
			def org_db_built = ORG_DB_BUILD().map({ d -> d.dir })
			org_db_ch = org_db_built
			built_ch = built_ch.mix(org_db_built.map({ d -> record(name: 'org_db', dir: d) }))
		}
	}

	// AMRFinderPlus
	def amrfinderplus_ch: Value<Path?> = channel.value(null)
	if (need_amrfinderplus) {
		if (amrfinderplus_dir != null) {
			amrfinderplus_ch = channel.value(amrfinderplus_dir)
		} else if (is_archive) {
			amrfinderplus_ch = root.map({ d -> (d as Path).resolve('amrfinderplus') })
		} else if (dir != null && dir.resolve('amrfinderplus').exists()) {
			amrfinderplus_ch = channel.value(dir.resolve('amrfinderplus'))
		} else {
			def amrfinderplus_built = AMRFINDERPLUS_UPDATE()
			amrfinderplus_ch = amrfinderplus_built
			built_ch = built_ch.mix(amrfinderplus_built.map({ d -> record(name: 'amrfinderplus', dir: d) }))
		}
	}

	emit:
	checkm2: Value<Path?>       = checkm2_ch          // CheckM2 database folder
	org_db: Value<Path?>        = org_db_ch           // orgfinder database folder, as built by ORG_DB_BUILD
	amrfinderplus: Value<Path?> = amrfinderplus_ch    // AMRFinderPlus database folder
	built: Channel<BuiltDb>     = built_ch            // the databases built by this run
}
