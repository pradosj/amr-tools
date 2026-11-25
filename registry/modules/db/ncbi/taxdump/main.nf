

process NCBI_TAXDUMP_DOWNLOAD {
	  output:
			path('taxdump/')
	  script:
	    """
	    mkdir -p taxdump \
	      && curl -kL https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz \
	      | tar -C ./taxdump -zxf -
	    """
}

