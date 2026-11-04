nextflow.enable.types = true

include { CHECKM2_DB           } from './modules/checkm2/db'
include { AMRFINDERPLUS_UPDATE } from './modules/amrfinderplus/update'
include { UNTAR                } from './modules/tar/untar.nf'
include { ORG_DB_BUILD         } from './subworkflows/org_db/main.nf'


// A database built by DATABASES, published under its name in the databases folder
record BuiltDb {
	name: String    // path relative to the databases folder (checkm2, org_db, amrfinderplus by default)
	dir: Path
}


// Gather the databases in a databases folder, as sub-folders checkm2/, org_db/ and amrfinderplus/ (links to the databases)
process DB_FOLDER {
	label 'local'

	input:
	checkm2: Path
	org_db: Path
	amrfinderplus: Path

	stage:
	stageAs(checkm2, 'checkm2')
	stageAs(org_db, 'org_db')
	stageAs(amrfinderplus, 'amrfinderplus')

	output:
	file('db')

	script:
	"""
	mkdir db
	mv checkm2 org_db amrfinderplus db/
	"""
}


// A database of the pipelines
enum Db {
	CHECKM2,
	ORG_DB,
	AMRFINDERPLUS
}


// Databases folder of the pipelines (params.db). Fields are nullable: defaults are used when null.
record DbParams {
	basedir: String?          // databases folder (default: results/db), or tar archive of it (path or URL)
	checkm2: String?          // CheckM2 database, relative to basedir (default: checkm2)
	org_db: String?           // orgfinder database, relative to basedir (default: org_db)
	amrfinderplus: String?    // AMRFinderPlus database, relative to basedir (default: amrfinderplus)
}


// Databases of the pipelines, from the databases folder `opts.basedir` (a database missing from it is built,
// and published by the pipelines under its relative path in the db/ output folder, i.e. results/db by default),
// or from a tar archive of it (path or URL), that must hold all the needed databases.
// A database that is not needed is null.
workflow DATABASES {
	take:
	opts: DbParams
	needs: Set<Db>    // the databases needed by the pipeline (the others are null)

	main:
	def basedir = opts.basedir ?: 'results/db'
	def checkm2_path = opts.checkm2 ?: 'checkm2'
	def org_db_path = opts.org_db ?: 'org_db'
	def amrfinderplus_path = opts.amrfinderplus ?: 'amrfinderplus'
	def is_archive = basedir ==~ /.+\.(tar|tar\.gz|tgz|tar\.bz2|tar\.xz)$/
	def dir = is_archive ? null : file(basedir)
	def root: Value<Path?> = channel.value(null)
	if (is_archive && !needs.isEmpty()) {
		root = UNTAR(channel.value(record(archive: file(basedir), dir_name: 'db')))
	}
	def built_ch: Channel<BuiltDb> = channel.empty()

	// CheckM2
	def checkm2_ch: Value<Path?> = channel.value(null)
	if (Db.CHECKM2 in needs) {
		if (is_archive) {
			checkm2_ch = root.map({ d -> (d as Path).resolve(checkm2_path) })
		} else if (dir.resolve(checkm2_path).exists()) {
			checkm2_ch = channel.value(dir.resolve(checkm2_path))
		} else {
			def checkm2_built = CHECKM2_DB(checkm2_path)
			checkm2_ch = checkm2_built
			built_ch = built_ch.mix(checkm2_built.map({ d -> record(name: checkm2_path, dir: d) }))
		}
	}

	// orgfinder
	def org_db_ch: Value<Path?> = channel.value(null)
	if (Db.ORG_DB in needs) {
		if (is_archive) {
			org_db_ch = root.map({ d -> (d as Path).resolve(org_db_path) })
		} else if (dir.resolve(org_db_path).exists()) {
			org_db_ch = channel.value(dir.resolve(org_db_path))
		} else {
			def org_db_built = ORG_DB_BUILD().map({ d -> d.dir })
			org_db_ch = org_db_built
			built_ch = built_ch.mix(org_db_built.map({ d -> record(name: org_db_path, dir: d) }))
		}
	}

	// AMRFinderPlus
	def amrfinderplus_ch: Value<Path?> = channel.value(null)
	if (Db.AMRFINDERPLUS in needs) {
		if (is_archive) {
			amrfinderplus_ch = root.map({ d -> (d as Path).resolve(amrfinderplus_path) })
		} else if (dir.resolve(amrfinderplus_path).exists()) {
			amrfinderplus_ch = channel.value(dir.resolve(amrfinderplus_path))
		} else {
			def amrfinderplus_built = AMRFINDERPLUS_UPDATE(amrfinderplus_path)
			amrfinderplus_ch = amrfinderplus_built
			built_ch = built_ch.mix(amrfinderplus_built.map({ d -> record(name: amrfinderplus_path, dir: d) }))
		}
	}

	emit:
	checkm2: Value<Path?>       = checkm2_ch          // CheckM2 database folder
	org_db: Value<Path?>        = org_db_ch           // orgfinder database folder, as built by ORG_DB_BUILD
	amrfinderplus: Value<Path?> = amrfinderplus_ch    // AMRFinderPlus database folder
	built: Channel<BuiltDb>     = built_ch            // the databases built by this run, named by their path relative to basedir
}
