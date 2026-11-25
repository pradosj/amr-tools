#!/usr/bin/awk -f

#Usage: fq_subsample.awk -v bp_limit=1000000 read_1.fq.gz read_2.fq.gz

BEGIN{
	bp=0
	
	# Read a line from each input file
	current_line = 1
	has_line = 1
	while(has_line) {

		# Read a line from each input FASTQ
		for(i=1;i<length(ARGV);i++){has_line = has_line && (("gzip -dc " ARGV[i]) | getline line[i "_" current_line%4])}
		
		if (!has_line){exit 0}
		
		# Every batch of 4 lines
		if (current_line % 4 == 0) {
		  for(i=1;i<length(ARGV);i++) {
		    if (line[i "_" 1]!~/^@/) {print "Expected @ at line " current_line " of file " ARGV[i] > "/dev/stderr";exit 1}
		    if (line[i "_" 3]!~/^[+]/) {print "Expected + at line " current_line " of file " ARGV[i] > "/dev/stderr";exit 1}
		    bp += length(line[i "_" 2])
		    ocmd = "gzip --fast > subsampled_reads_" i ".fastq.gz"
		    print line[i "_" 1] | ocmd
		    print line[i "_" 2] | ocmd
		    print line[i "_" 3] | ocmd
		    print line[i "_" 0] | ocmd
		  }
		  if (bp>=bp_limit){exit 0}
		}
		current_line += 1
	}
}


