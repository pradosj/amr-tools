

include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

process GTF2GFF {
	  container "registry.gitlab.unige.ch/amr-genomics/rscript:v3"
    memory '4 GB'
    cpus 1
    time '30 min'
    input:
    		tuple(val(meta),path("ref.gtf"))
    		path(extra)
    output:
        tuple(val(meta),path('ref.gff'))
    script:
				"""
				#!/usr/bin/env Rscript
				source("assets/lib_gtf2gff.R")
				rtracklayer::import.gff2("ref.gtf") %>% 
					gtf2gff() %>% 
					rtracklayer::export.gff3("ref.gff")
				"""
}


workflow {
	main:
		// Validate parameters and print summary of supplied ones
		validateParameters()
		log.info(paramsSummaryLog(workflow))

		def gtf_ch = channel.fromPath(params.gtf).map({f ->
				def id = f.name.replaceAll(/\.(gtf)$/,'')
				[[sample_id:id],f]
		})
		
		GTF2GFF(gtf_ch,file("${projectDir}/assets"))

	publish:
		gff = GTF2GFF.out
}

output {
	gff {
		path { m,x -> x >> "${m.sample_id}.gff"}
	}
}

