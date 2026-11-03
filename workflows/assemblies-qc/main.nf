nextflow.enable.types = true

include { MINIMAP2_ALIGN_ONT } from './modules/minimap2/align_ont'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_LONG  } from './modules/samtools/stats'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_SHORT } from './modules/samtools/stats'
include { BWA_MEM            } from './modules/bwa/mem'
include { BWA_INDEX          } from './modules/bwa/index'
include { CHECKM2_DB         } from './modules/checkm2/db'
include { CHECKM2_PREDICT    } from './modules/checkm2/predict'


// Options of CheckM2 (params.assemblies_qc.checkm2)
record Checkm2Params {
	skip: Boolean?
	db: Path?                  // CheckM2 database folder, downloaded when null
	args: String?
}

// Options of ASSEMBLIES_QC (params.assemblies_qc)
record AssembliesQcParams {
	checkm2: Checkm2Params?
}

// Assembly to check. meta holds both sample_id and assembler_name.
record Genome {
	meta: Map
	fasta: Path
	gff: Path?
	info: Path?
}

// Illumina reads of the sample of an assembly (same meta as the assembly): read2 is null for single-end reads
record ShortReads {
	meta: Map?
	read1: Path
	read2: Path?
}

// ONT reads of the sample of an assembly (same meta as the assembly): one FASTQ
record LongReads {
	meta: Map?
	reads: Path
}


// Gather the QC statistics of an assembly into an RDS file
process COMBINE_QC_STATS {
	container "registry.gitlab.unige.ch/amr-genomics/rscript:main"
	label 'rscript'
	memory 8.GB
	cpus 2
	time 30.m
	input:
	record(meta: Map, fasta: Path, long_stats: Path?, short_stats: Path?, checkm2: Path?)
	lib: Path    // lib_assembly_stats.R
	stage:
	stageAs(fasta, 'stats/assembly.fasta')
	stageAs(long_stats, 'stats/long_reads.bam.stats')
	stageAs(short_stats, 'stats/short_reads.bam.stats')
	stageAs(checkm2, 'stats/checkm2_quality_report.tsv')
	output:
	record(meta: meta, rds: file('qc.rds'))
	script:
	// meta is reported as [[sample_id], [assembler_name]]
	def meta_json = new groovy.json.JsonBuilder([meta.subMap(['sample_id']), meta.subMap(['assembler_name'])]).toString()
	"""
	#!/usr/bin/env Rscript
	source("${lib}")
	stats <- read_assembly_stats("./stats")
	stats[["meta"]] <- '${meta_json}'
	stats[["sample_id"]] <- '${meta.sample_id}'
	stats[["assembler_name"]] <- '${meta.assembler_name}'
	saveRDS(stats,"qc.rds")
	"""
	stub:
	"""
	touch qc.rds
	"""
}

// Aggregate the QC statistics of all the assemblies into a table
process ASSEMBLIES_MULTIQC {
	container "registry.gitlab.unige.ch/amr-genomics/rscript:main"
	label 'rscript'
	memory 8.GB
	cpus 2
	time 30.m
	input:
	rds: List<Path>
	lib: Path    // lib_assembly_stats.R
	stage:
	stageAs(rds, 'stats/qc*.rds')
	output:
	record(
		xlsx: file('assemblies_multiqc.xlsx'),
		txt: file('assemblies_multiqc.txt')
	)
	script:
	"""
	#!/usr/bin/env Rscript
	source("${lib}")
	A <- rlang::inject(aggregate_assembly_stats(!!!(c(${rds.collect({ f -> "'${f}'" }).join(",")}) |> map(readRDS) |> unname())))
	A |> write.table(file="assemblies_multiqc.txt",sep="\\t",quote=FALSE,row.names=FALSE)
	A |> openxlsx::write.xlsx("assemblies_multiqc.xlsx")
	"""
	stub:
	"""
	touch assemblies_multiqc.xlsx assemblies_multiqc.txt
	"""
}


workflow ASSEMBLIES_QC {
	take:
	opts: AssembliesQcParams
	assemblies: Channel<Genome>
	short_reads: Channel<ShortReads>
	long_reads: Channel<LongReads>

	main:
	def lib = file("${moduleDir}/assets/lib_assembly_stats.R")

	// Short reads alignment and statistics
	def short_bam_ch = BWA_MEM(BWA_INDEX(assemblies).join(short_reads, by: 'meta'))
	def short_stats_ch = SAMTOOLS_STATS_SHORT(short_bam_ch)

	// Long reads alignment and statistics
	def long_bam_ch = MINIMAP2_ALIGN_ONT(assemblies.join(long_reads, by: 'meta'))
	def long_stats_ch = SAMTOOLS_STATS_LONG(long_bam_ch)

	// Genome completeness
	def checkm2_db_ch = channel.empty()
	def checkm2_ch = channel.empty()
	if (!opts.checkm2?.skip) {
		def db = opts.checkm2?.db ? channel.value(opts.checkm2.db) : CHECKM2_DB()
		checkm2_db_ch = channel.empty().mix(db)
		checkm2_ch = CHECKM2_PREDICT(assemblies.combine(db: db).map({ r -> r + record(args: opts.checkm2?.args) }))
	}

	//TODO: RUN VCF_LONG
	//TODO: RUN VCF_SHORT
	//TODO: HTML_REPORT()
	//TODO: CHARACTERIZE_UNMAPPED_READS
	def qc_ch = COMBINE_QC_STATS(
		assemblies
			.join(long_stats_ch.map({ r -> record(meta: r.meta, long_stats: r.stats) }), by: 'meta', remainder: true)
			.join(short_stats_ch.map({ r -> record(meta: r.meta, short_stats: r.stats) }), by: 'meta', remainder: true)
			.join(checkm2_ch.map({ r -> record(meta: r.meta, checkm2: r.report) }), by: 'meta', remainder: true),
		lib
	)

	// no table when there is no assembly
	def multiqc_ch = ASSEMBLIES_MULTIQC(
		qc_ch.map({ r -> r.rds }).collect().flatMap({ fs -> fs.isEmpty() ? [] : [fs.toList()] }),
		lib
	)

	emit:
	long_alignments: Channel<Record>  = long_bam_ch       // records with fields meta, bam, bai
	short_alignments: Channel<Record> = short_bam_ch      // records with fields meta, bam, bai
	assembly_qc: Channel<Record>      = qc_ch             // records with fields meta, rds
	multiqc: Channel<Record>          = multiqc_ch        // records with fields txt, xlsx
	checkm2: Channel<Record>          = checkm2_ch        // records with fields meta, report
	checkm2_db: Channel<Path>         = checkm2_db_ch
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
