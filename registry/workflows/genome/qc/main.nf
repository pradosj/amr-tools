nextflow.enable.types = true

include { MINIMAP2_ALIGN_BAM } from './modules/minimap2/align_bam'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_LONG  } from './modules/samtools/stats'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_SHORT } from './modules/samtools/stats'
include { SAMTOOLS_COVERAGE as SAMTOOLS_COVERAGE_LONG  } from './modules/samtools/coverage'
include { SAMTOOLS_COVERAGE as SAMTOOLS_COVERAGE_SHORT } from './modules/samtools/coverage'
include { BWA_MEM            } from './modules/bwa/mem'
include { BWA_INDEX          } from './modules/bwa/index'
include { CHECKM2_PREDICT    } from './modules/checkm2/predict'
include { MetaPath           } from './subworkflows/input/main.nf'


// Options of CheckM2 (params.assemblies_qc.checkm2)
record Checkm2Params {
	skip: Boolean?
	args: String?
}

// Options of GENOME_QC (params.assemblies_qc)
record AssembliesQcParams {
	checkm2: Checkm2Params?
}

// Illumina reads of the sample of an assembly (same meta as the assembly): read2 is null for single-end reads
record ShortReads {
	meta: Record?
	read1: Path
	read2: Path?
}

// ONT reads of the sample of an assembly (same meta as the assembly): one FASTQ
record LongReads {
	meta: Record?
	reads: Path
}


// Gather the QC statistics of an assembly into an RDS file
process COMBINE_QC_STATS {
	container "registry.gitlab.unige.ch/amr-genomics/rscript:v3"
	label 'rscript'
	memory 8.GB
	cpus 2
	time 30.m
	input:
	record(meta: Record, fasta: Path, long_stats: Path?, short_stats: Path?, long_cov: Path?, short_cov: Path?, checkm2: Path?)
	assets: Path    // the assets folder of the workflow (lib_io.R, lib_assembly_stats.R)
	stage:
	stageAs(assets, 'assets')
	stageAs(fasta, 'stats/assembly.fasta')
	stageAs(long_stats, 'stats/long_reads.bam.stats')
	stageAs(short_stats, 'stats/short_reads.bam.stats')
	stageAs(long_cov, 'stats/long_reads.coverage.tsv')
	stageAs(short_cov, 'stats/short_reads.coverage.tsv')
	stageAs(checkm2, 'stats/checkm2_quality_report.tsv')
	output:
	record(meta: meta, rds: file('qc.rds'))
	script:
	// meta is reported as [[sample_id], [assembler_name]]
	def meta_json = new groovy.json.JsonBuilder([[sample_id: meta.id], [assembler_name: meta.assembler_name]]).toString()
	"""
	#!/usr/bin/env Rscript
	source("assets/lib_io.R")
	source("assets/lib_assembly_stats.R")
	stats <- read_assembly_stats("./stats")
	stats[["meta"]] <- '${meta_json}'
	stats[["sample_id"]] <- '${meta.id}'
	stats[["assembler_name"]] <- '${meta.assembler_name}'
	saveRDS(stats,"qc.rds")
	"""
	stub:
	"""
	touch qc.rds
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
	rds: List<Path>
	assets: Path    // the assets folder of the workflow (lib_io.R, lib_assembly_stats.R)
	stage:
	stageAs(assets, 'assets')
	stageAs(rds, 'stats/qc*.rds')
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
	A <- rlang::inject(aggregate_assembly_stats(!!!(c(${rds.collect({ f -> "'${f}'" }).join(",")}) |> map(readRDS) |> unname())))
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
	opts: AssembliesQcParams
	checkm2_db: Value<Path?>          // CheckM2 database folder (null when checkm2 is skipped)
	assemblies: Channel<MetaPath>     // the FASTA of each assembly, meta holding id (sample id) and assembler_name
	short_reads: Channel<ShortReads>
	long_reads: Channel<LongReads>

	main:
	def assets = file("${moduleDir}/assets")

	// Short reads alignment on the BWA index of the assembly, and statistics
	def index_ch = BWA_INDEX(assemblies)
		.map({ i -> record(meta: i.meta, ref: i) })
	def short_bam_ch = BWA_MEM(short_reads.join(index_ch, by: 'meta'))
	def short_stats_ch = SAMTOOLS_STATS_SHORT(short_bam_ch.bam)
	def short_cov_ch = SAMTOOLS_COVERAGE_SHORT(short_bam_ch.bam)

	// Long reads alignment on the assembly, and statistics
	def long_bam_ch = MINIMAP2_ALIGN_BAM(
		long_reads
			.map({ r -> record(meta: r.meta, path: r.reads) })
			.join(assemblies.map({ r -> record(meta: r.meta, ref: r, args: '-x map-ont') }), by: 'meta')
	)
	def long_stats_ch = SAMTOOLS_STATS_LONG(long_bam_ch.bam)
	def long_cov_ch = SAMTOOLS_COVERAGE_LONG(long_bam_ch.bam)

	// Genome completeness
	def checkm2_ch = channel.empty()
	if (!opts.checkm2?.skip) {
		checkm2_ch = CHECKM2_PREDICT(assemblies.combine(db: checkm2_db, args: opts.checkm2?.args))
	}

	//TODO: RUN VCF_LONG
	//TODO: RUN VCF_SHORT
	//TODO: HTML_REPORT()
	//TODO: CHARACTERIZE_UNMAPPED_READS
	def qc_ch = COMBINE_QC_STATS(
		assemblies
			.map({ r -> record(meta: r.meta, fasta: r.path) })
			.join(long_stats_ch.map({ r -> record(meta: r.meta, long_stats: r.path) }), by: 'meta', remainder: true)
			.join(short_stats_ch.map({ r -> record(meta: r.meta, short_stats: r.path) }), by: 'meta', remainder: true)
			.join(long_cov_ch.map({ r -> record(meta: r.meta, long_cov: r.path) }), by: 'meta', remainder: true)
			.join(short_cov_ch.map({ r -> record(meta: r.meta, short_cov: r.path) }), by: 'meta', remainder: true)
			.join(checkm2_ch.map({ r -> record(meta: r.meta, checkm2: r.path) }), by: 'meta', remainder: true),
		assets
	)

	// no table when there is no assembly
	def multiqc_ch = GENOME_MULTIQC(
		qc_ch.map({ r -> r.rds }).collect().flatMap({ fs -> fs.isEmpty() ? [] : [fs.toList()] }),
		assets
	)

	emit:
	long_alignments: Channel<Record>  = long_bam_ch.outputs     // records with fields meta, meta_ref, bam, bai
	short_alignments: Channel<Record> = short_bam_ch.outputs    // records with fields meta, meta_ref, bam, bai
	assembly_qc: Channel<Record>      = qc_ch             // records with fields meta, rds
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
