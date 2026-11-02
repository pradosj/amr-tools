
include { MINIMAP2_ALIGN_ONT } from './modules/minimap2/align_ont'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_LONG  } from './modules/samtools/stats'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_SHORT } from './modules/samtools/stats'
include { BWA_MEM            } from './modules/bwa/mem'
include { BWA_INDEX          } from './modules/bwa/index'
include { ORGANIZE_FILES     } from './modules/organize_files'
include { CHECKM2_DB         } from './modules/checkm2/db'
include { CHECKM2_PREDICT    } from './modules/checkm2/predict'



process COMBINE_QC_STATS {
  container "registry.gitlab.unige.ch/amr-genomics/rscript:main"
  memory '8 GB'
  cpus 2
  time '30 min'
  input:
		tuple val(meta), path("stats")
		path("assets") 
  output:
    tuple val(meta), path('qc.rds'), emit: 'rds'
  script:
		def builder = new groovy.json.JsonBuilder(meta)
	  """
		#!/usr/bin/env Rscript
		source("assets/lib_assembly_stats.R")
		stats <- read_assembly_stats("./stats")
		stats[["meta"]] <- '${builder.toString()}'
		stats[["sample_id"]] <- '${meta[0].sample_id}'
		stats[["assembler_name"]] <- '${meta[1].assembler_name}'
		saveRDS(stats,"qc.rds")
		"""
}

process ASSEMBLIES_MULTIQC {
  container "registry.gitlab.unige.ch/amr-genomics/rscript:main"
  memory '8 GB'
  cpus 2
  time '30 min'
  input:
		path(rds,name: 'stats/qc*.rds')
		path("assets") 
  output:
    path('assemblies_multiqc.xlsx'), emit: xlsx
    path('assemblies_multiqc.txt'),  emit: txt
  script:
		"""
		#!/usr/bin/env Rscript
		source("assets/lib_assembly_stats.R")
		A <- rlang::inject(aggregate_assembly_stats(!!!(c(${[rds].flatten().collect({f -> "'" + f + "'"}).join(",")}) |> map(readRDS) |> unname()))) 
		A |> write.table(file="assemblies_multiqc.txt",sep="\\t",quote=FALSE,row.names=FALSE)
		A |> openxlsx::write.xlsx("assemblies_multiqc.xlsx")
		"""
}

workflow ASSEMBLIES_QC {
	take:
		opts
		fa_ch
		fqs_ch
		fql_ch
	main:
	
		// Short reads alignment and statistics
		BWA_MEM(BWA_INDEX(fa_ch).join(fqs_ch))
		SAMTOOLS_STATS_SHORT(BWA_MEM.out.bam)
		
		// Long reads alignment and statistics
		MINIMAP2_ALIGN_ONT(fa_ch.join(fql_ch))
		SAMTOOLS_STATS_LONG(MINIMAP2_ALIGN_ONT.out.bam)
		
		// Genome completness
		def checkm2_db_ch = channel.empty()
		def checkm2_ch = channel.empty()
		if (!opts.checkm2.skip) {
			if (opts.checkm2.db) {
				checkm2_db_ch = channel.FromPath(opts.checkm2.db)
			} else {
				checkm2_db_ch = CHECKM2_DB()
			}
			CHECKM2_PREDICT(fa_ch,checkm2_db_ch)
			checkm2_ch = CHECKM2_PREDICT.out
		}
		
		//TODO: RUN VCF_LONG
		//TODO: RUN VCF_SHORT
		//TODO: HTML_REPORT()
		//TODO: CHARACTERIZE_UNMAPPED_READS
		fa_ch
			.join(SAMTOOLS_STATS_LONG.out,remainder:true)
			.join(SAMTOOLS_STATS_SHORT.out,remainder:true)
			.join(checkm2_ch,remainder:true)
			.map({meta,x1,x2,x3,x4 -> [meta,[[x1,"assembly.fasta"],[x2,"long_reads.bam.stats"],[x3,"short_reads.bam.stats"],[x4,"checkm2_quality_report.tsv"]].findAll({x,_y -> x})]})
			| ORGANIZE_FILES
		
		COMBINE_QC_STATS(ORGANIZE_FILES.out,"${moduleDir}/assets")

		def assembly_multiqc_txt_ch = channel.empty()
		def assembly_multiqc_xlsx_ch = channel.empty()
		def rds = COMBINE_QC_STATS.out.rds.map({_m,x -> x}).toList()
		if (rds) {
			ASSEMBLIES_MULTIQC(rds,"${moduleDir}/assets")
			assembly_multiqc_txt_ch  = ASSEMBLIES_MULTIQC.out.txt
			assembly_multiqc_xlsx_ch = ASSEMBLIES_MULTIQC.out.xlsx
		}

	emit:
		long_bam        = MINIMAP2_ALIGN_ONT.out.bam
		long_bai        = MINIMAP2_ALIGN_ONT.out.bai
		long_bam_stats  = SAMTOOLS_STATS_LONG.out
		long_vcf        = channel.empty()
		
		short_bam       = BWA_MEM.out.bam
		short_bai       = BWA_MEM.out.bai
		short_bam_stats = SAMTOOLS_STATS_SHORT.out
		short_vcf       = channel.empty()
		
		assembly_qc           = COMBINE_QC_STATS.out.rds
		assembly_multiqc_txt  = assembly_multiqc_txt_ch
		assembly_multiqc_xlsx = assembly_multiqc_xlsx_ch
		
		checkm2    = checkm2_ch
		checkm2_db = checkm2_db_ch
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