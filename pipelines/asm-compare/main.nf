

include { MINIMAP2_ALIGN_ASM5   } from './modules/minimap2/align_asm5'
include { BCFTOOLS_ASM5_MPILEUP } from './modules/bcftools/asm5_mpileup'
include { QUARTO_RENDER         } from './modules/quarto/render'
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'


process ROTATE_FASTA {
	  container "registry.gitlab.unige.ch/amr-genomics/rscript:v2"
    memory '8 GB'
    cpus 1
    time '30 min'
    input:
    		tuple(val(meta),path("query.fasta"),path("aln.bam"))
    		path(extra)
    output:
        tuple(val(meta),path('rotated_query.fasta'))
    script:
				"""
				#!/usr/bin/env Rscript
				source("assets/lib_bam.R")
				rotate_queries("query.fasta","aln.bam") |>
					writeXStringSet("rotated_query.fasta")
				"""
}

workflow ROTATE_QUERIES {
	take:
		query_ch  // [meta, path('ref.fasta'), path('qry.fasta')]
	main:
		query_ch | MINIMAP2_ALIGN_ASM5
		
		def rotated_query_ch = ROTATE_FASTA(
			query_ch
			.map({m,_r,q -> [m,q]})
			.combine(MINIMAP2_ALIGN_ASM5.out.bam,by:0),
			file("${projectDir}/assets")
		)
		.combine(query_ch,by:0)
		.map({m,rq,r,_q -> [m,r,rq]})
	emit:
		rotated_queries = rotated_query_ch
}


workflow {
	main:
		// Validate parameters and print summary of supplied ones
		validateParameters()
		log.info(paramsSummaryLog(workflow))

		def ref = file(params.ref)
		def ref_gff = params.ref_gff ? file(params.ref_gff) : []
		def query_ch = channel.fromPath(params.query).map({f ->
				def id = f.name.replaceAll(/\.(fasta|fna|fa)$/,'')
				[[sample_id:id],ref,f]
		})

		if (!params.no_query_rotation) {
			query_ch = ROTATE_QUERIES(query_ch)
		}

		// Map (rotated) queries to REF
		query_ch
		| MINIMAP2_ALIGN_ASM5
		
		MINIMAP2_ALIGN_ASM5.out.bam
		| map({meta,bam -> [meta,ref,bam,ref_gff]})
		| BCFTOOLS_ASM5_MPILEUP

		QUARTO_RENDER(
			MINIMAP2_ALIGN_ASM5.out.bam
				.combine(BCFTOOLS_ASM5_MPILEUP.out.txt,by:0)
				.map({meta,bam,txt -> [meta,[bam,txt],"-P 'bam_file:${bam}' -P 'txt_file:${txt}'"]}),
			file("${projectDir}/assets/bam_report.qmd"),
			file("${projectDir}/assets")
		)

	publish:
		rotated_query_fa = query_ch.filter({!params.no_query_rotation}).map({m,_r,q -> [m,q]})
		bam = MINIMAP2_ALIGN_ASM5.out.bam
		bai = MINIMAP2_ALIGN_ASM5.out.bai
		mut_vcf = BCFTOOLS_ASM5_MPILEUP.out.vcf
		mut_vcf_csi = BCFTOOLS_ASM5_MPILEUP.out.vcf_csi
		mut_txt = BCFTOOLS_ASM5_MPILEUP.out.txt
		html = QUARTO_RENDER.out.html
}

output {
	rotated_query_fa {
		path { m,x -> x >> "${m.sample_id}.rotated.fasta"}
	}
	bam {
		path { m,x -> x >> "${m.sample_id}.bam"}
	}
	bai {
		path { m,x -> x >> "${m.sample_id}.bai"}
	}
	mut_vcf {
		path { m,x -> x >> "${m.sample_id}.vcf.gz"}
	}
	mut_vcf_csi {
		path { m,x -> x >> "${m.sample_id}.vcf.gz.csi"}
	}
	mut_txt {
		path { m,x -> x >> "${m.sample_id}.txt"}
	}	
	html {
		path {  m,x -> x >> "${m.sample_id}.html"}
	}
	
}




