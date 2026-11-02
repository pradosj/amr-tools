
process ORGANIZE_FILES {
    memory '2 GB'
    cpus 1
    time '1h'
    executor 'local'
    input:
    	tuple(val(meta),val(file_pairs)) // liste de [[path(src), target]]
    output:
    	tuple(val(meta),path("output",type:'dir'))
    script:
    	def cmd = file_pairs
    		.findAll({src,_target -> src})
    		.collect({src,target ->
					def dest_path = "output/${target}"
					if (src instanceof java.nio.file.Path) {	
					    return "mkdir -p \$(dirname '${dest_path}') && ln -s '${src}' '${dest_path}'".stripIndent()
					} else {
							return "mkdir -p \$(dirname '${dest_path}') && echo '${src}' > '${dest_path}'".stripIndent()
					}
	    	})
	    """
	    trap '' PIPE
	    mkdir -p output
	    ${cmd.join('\n')}
	    """
}


