
process RSCRIPT {
	  container "registry.gitlab.unige.ch/amr-genomics/rscript:main"
	  label 'rscript'
    memory '8 GB'
    cpus 1
    time '30 min'
    input:
    		tuple(val(meta), path(inputs))
    		each path('rscript.R')
    		each path('extra/*')
    output:
        tuple(val(meta), path('output'))
    script:
				"""
				Rscript rscript.R
				"""
}
