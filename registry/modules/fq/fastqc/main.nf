process FASTQC {
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    cpus 3
    memory '5 GB'
    time '1h'
    input:
	    tuple val(meta), path(reads), val(id)
    output:
	    tuple val(meta), path("*_fastqc.html"), emit: html
	    tuple val(meta), path("*_fastqc.zip"), emit: zip
    script:
	    """
	    gzip -dc ${reads} | gzip > ${id}.fastq.gz
	    fastqc \\
	        ${task.ext.args?:''} \\
	        --threads ${task.cpus} \\
	        --memory 5000 \\
	        ${id}.fastq.gz
	    rm ${id}.fastq.gz
	    """
	  stub:
	  	"""
	  	touch ${id}_fastqc.html
	  	touch ${id}_fastqc.zip
	  	"""
}
