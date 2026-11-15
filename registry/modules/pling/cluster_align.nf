process PLING_CLUSTER_ALIGN {
    container 'docker.io/staphb/pling:3.0.1'
    memory '10 GB'
    cpus 8
    time '1h'
    input:
	    tuple val(meta), path('plasmids.fasta')
    output:
	    tuple val(meta), path('pling_out/',type:'dir')
    script:
	    """
	    # Break the given multi-fasta into multiple fasta files
	    mkdir -p plasmids
      awk '/^>/{if(f) close(f); id=id+1; f="plasmids/" id ".fasta"} {print > f}' plasmids.fasta
      
      # Run pling
      ls plasmids/*.fasta > input.txt
	    pling cluster align input.txt pling_out ${task.ext.args?:''} --cores ${task.cpus}
	    """
		stub:
			"""
			mkdir -p pling_out
			"""
}

