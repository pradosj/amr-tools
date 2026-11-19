nextflow.enable.types = true

include { MINIMAP2_ALIGN_BAM } from './modules/minimap2/align_bam'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_LONG  } from './modules/samtools/stats'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_SHORT } from './modules/samtools/stats'
include { SAMTOOLS_COVERAGE as SAMTOOLS_COVERAGE_LONG  } from './modules/samtools/coverage'
include { SAMTOOLS_COVERAGE as SAMTOOLS_COVERAGE_SHORT } from './modules/samtools/coverage'
include { BWA_MEM            } from './modules/bwa/mem'
include { BWA_INDEX          } from './modules/bwa/index'
include { CHECKM2_PREDICT    } from './modules/checkm2/predict'
include { QUARTO_RENDER      } from './modules/quarto/render'
include { MetaPath; ShortReads; LongReads } from './subworkflows/input/main.nf'


// Options of CheckM2 (params.genome_qc.checkm2 of the pipelines)
record Checkm2Params {
	skip: Boolean?
	args: String?
	db: Path?        // CheckM2 database folder
}

// Options of the QC of the assemblies (params.genome_qc of the pipelines)
record GenomeQcParams {
	checkm2: Checkm2Params?
}

// An assembly and the reads of its sample (meta holding id (sample id) and assembler_name)
record GenomeQcInput {
	meta: Record?
	path: Path                  // FASTA of the assembly
	short_reads: ShortReads?    // Illumina reads (their meta is not used), null when there are none
	long_reads: LongReads?      // ONT reads (their meta is not used), null when there are none
}

// Gather the QC statistics of an assembly into a JSON file (read_assembly_stats() of lib_assembly_stats.R)
process MERGE_QC_STATS {
	container "registry.gitlab.unige.ch/amr-genomics/rscript:v3"
	label 'rscript'
	memory 8.GB
	cpus 2
	time 30.m
	input:
	record(meta: Record?, path: Path, long_stats: Path?, short_stats: Path?, long_cov: Path?, short_cov: Path?, checkm2: Path?)
	assets: Path    // the assets folder of the workflow (lib_io.R, lib_assembly_stats.R)
	stage:
	stageAs(assets, 'assets')
	stageAs(path, 'stats/assembly.fasta')
	stageAs(long_stats, 'stats/long_reads.bam.stats')
	stageAs(short_stats, 'stats/short_reads.bam.stats')
	stageAs(long_cov, 'stats/long_reads.coverage.tsv')
	stageAs(short_cov, 'stats/short_reads.coverage.tsv')
	stageAs(checkm2, 'stats/checkm2_quality_report.tsv')
	output:
	record(meta: meta, json: file('qc.json'))
	script:
	"""
	#!/usr/bin/env Rscript
	source("assets/lib_io.R")
	source("assets/lib_assembly_stats.R")
	qc <- merge_genome_qc_stats("./stats")
	qc[['meta']][['id']] <- '${meta.id}'
	jsonlite::write_json(qc, "qc.json", auto_unbox = TRUE, na = "null", null = "null", digits = NA, pretty = TRUE)
	"""
	stub:
	"""
	touch qc.json
	"""
}

// Aggregate the QC statistics of all the assemblies into a table
process GENOME_MULTIQC {
	container "registry.gitlab.unige.ch/amr-genomics/rscript:v3"
	label 'rscript'
	memory 8.GB
	cpus 2
	time 30.m
	input:
	json: List<Path>    // the qc.json of each assembly
	assets: Path    // the assets folder of the workflow (lib_io.R, lib_assembly_stats.R)
	stage:
	stageAs(assets, 'assets')
	stageAs(json, 'stats/qc*.json')
	output:
	record(
		xlsx: file('assemblies_multiqc.xlsx'),
		txt: file('assemblies_multiqc.txt')
	)
	script:
	"""
	#!/usr/bin/env Rscript
	source("assets/lib_io.R")
	source("assets/lib_assembly_stats.R")
	A <- aggregate_assembly_stats(fs::dir_ls("stats", glob = "*.json"))
	A |> write.table(file="assemblies_multiqc.txt",sep="\\t",quote=FALSE,row.names=FALSE)
	A |> openxlsx::write.xlsx("assemblies_multiqc.xlsx")
	"""
	stub:
	"""
	touch assemblies_multiqc.xlsx assemblies_multiqc.txt
	"""
}


workflow GENOME_QC {
	take:
	genomes: Channel<GenomeQcInput> // each assembly with the reads of its sample
	opts: Value<GenomeQcParams?>    // options, with the CheckM2 database (CheckM2 skipped when checkm2.skip is set or its db is null)

	main:
	def assets = file("${moduleDir}/assets")

	// Short reads alignment on the BWA index of the assembly (indexed only when there are short reads), and statistics
	def with_short_ch = genomes.filter({ g -> g.short_reads != null })
	def index_ch = BWA_INDEX(with_short_ch.map({ g -> record(meta: g.meta, path: g.path) }))
		.map({ i -> record(meta: i.meta, ref: i) })
	def short_bam_ch = BWA_MEM(
		with_short_ch
			.map({ g -> record(meta: g.meta, read1: g.short_reads.read1, read2: g.short_reads.read2, args: null) })
			.join(index_ch, by: 'meta')
	)
	def short_stats_ch = SAMTOOLS_STATS_SHORT(short_bam_ch.bam)
	def short_cov_ch = SAMTOOLS_COVERAGE_SHORT(short_bam_ch.bam)

	// Long reads alignment on the assembly, and statistics
	def long_bam_ch = MINIMAP2_ALIGN_BAM(
		genomes
			.filter({ g -> g.long_reads != null })
			.map({ g -> record(meta: g.meta, path: g.long_reads.reads, ref: record(meta: g.meta, path: g.path), args: '-x map-ont') })
	)
	def long_stats_ch = SAMTOOLS_STATS_LONG(long_bam_ch.bam)
	def long_cov_ch = SAMTOOLS_COVERAGE_LONG(long_bam_ch.bam)

	// Genome completeness
	def checkm2_ch = CHECKM2_PREDICT(
		genomes
			.combine(opts: opts)
			.filter({ r -> !r.opts?.checkm2?.skip && r.opts?.checkm2?.db != null })
			.map({ r -> record(meta: r.meta, path: r.path, db: r.opts.checkm2.db, args: r.opts.checkm2.args) })
	)

	//TODO: RUN VCF_LONG
	//TODO: RUN VCF_SHORT
	//TODO: HTML_REPORT()
	//TODO: CHARACTERIZE_UNMAPPED_READS
	def qc_ch = MERGE_QC_STATS(
		genomes
			.join(long_stats_ch.map({ r -> record(meta: r.meta, long_stats: r.path) }), by: 'meta', remainder: true)
			.join(short_stats_ch.map({ r -> record(meta: r.meta, short_stats: r.path) }), by: 'meta', remainder: true)
			.join(long_cov_ch.map({ r -> record(meta: r.meta, long_cov: r.path) }), by: 'meta', remainder: true)
			.join(short_cov_ch.map({ r -> record(meta: r.meta, short_cov: r.path) }), by: 'meta', remainder: true)
			.join(checkm2_ch.map({ r -> record(meta: r.meta, checkm2: r.path) }), by: 'meta', remainder: true),
		assets
	)

	// HTML dashboard of each assembly
	def report_ch = QUARTO_RENDER(qc_ch.map({ r -> record(
		meta: r.meta,
		template: file("${moduleDir}/assets/genome_qc_report.qmd"),
		extra: assets,
		files: [r.json],
		args: "-P 'qc:${r.json.name}' -M 'title:Genome QC ${r.meta.id}__${r.meta.assembler_name}'"
	) }))

	// no table when there is no assembly
	def multiqc_ch = GENOME_MULTIQC(
		qc_ch.map({ r -> r.json }).collect().flatMap({ fs -> fs.isEmpty() ? [] : [fs.toList()] }),
		assets
	)

	emit:
	long_alignments: Channel<Record>  = long_bam_ch.outputs     // records with fields meta, meta_ref, bam, bai
	short_alignments: Channel<Record> = short_bam_ch.outputs    // records with fields meta, meta_ref, bam, bai
	assembly_qc: Channel<Record>      = qc_ch             // records with fields meta, json (qc.json)
	report: Channel<Record>           = report_ch         // records with fields meta, html (dashboard of each assembly)
	multiqc: Channel<Record>          = multiqc_ch        // records with fields txt, xlsx
	checkm2: Channel<MetaPath>        = checkm2_ch        // records with fields meta, path (quality_report.tsv)
}


/* Assembly stats
Short Reads
  - % properly paired
  - Number chimeric pair
  - Number of identified mutation in the VCF
Long Reads
  - % mapped
  - Number of identified mutation in the VCF
*/

/* Contigs stats
Short Reads
  - % mapped
  - % properly paired
  - Avg coverage
Long Reads
  - % mapped
  - Avg coverage
*/

/* Report
Short Reads
   - inter contigs links
*/


/*
process IGV_SCRIPT {
	input:
	  val(meta)
  output:
    tuple val(meta), path("load_in_igv.sh")
	script:
	"""
	file("${moduleDir}/assets/load_in_igv.sh")
	"""
}
*/
