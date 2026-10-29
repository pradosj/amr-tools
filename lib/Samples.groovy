
import nextflow.Channel


class Samples {
	
	static def fromCSV(csv_path,samplesheetToList) {
		return Channel.fromList(samplesheetToList(csv_path,"assets/schema/samples_sheets_schema.json"))
			.map({it[0]})
	}

	static def fromChannels(sr_ch,lr_ch,fa_ch) {
		def meta_ch = lr_ch.map({it[0]}).mix(sr_ch.map({it[0]}),fa_ch.map({it[0]})).unique()
		
		return meta_ch.ifEmpty(['__sentinel__',null])
			.join(lr_ch.ifEmpty(['__sentinel__',null]),failOnDuplicate:true,remainder:true)
			.join(sr_ch.ifEmpty(['__sentinel__',null]),failOnDuplicate:true,remainder:true)
			.join(fa_ch.ifEmpty(['__sentinel__',null]),failOnDuplicate:true,remainder:true)
			.filter({it[0]!='__sentinel__'})
			.map({it[0] + [
				reads_long:it[1],
				reads_short_1:it[2]?it[2][0]:null,
				reads_short_2:it[2]?it[2][1]:null,
				assembly_fasta: it[3]
			]})
	}

	static def fromParams(opts,samplesheetToList) {
		if (opts.csv) {
			return fromCSV(opts.csv,samplesheetToList)
		} else {
			def sr_ch = Channel.empty()
			def lr_ch = Channel.empty()
			def fa_ch = Channel.empty()
	  	if (opts.reads.illumina) {
	  		lr_ch = Channel.fromPath(opts.reads.illumina).map({
	  				def id = it.name.replaceAll(/\.(fastq\.gz|fq\.gz|bam|cram)$/,'')
	  				[[sample_id:id],it]
	  		})
	  	}
	  	if (opts.assemblies.fasta) {
				fa_ch = Channel.fromPath(opts.assemblies.fasta).map({
	  				def id = it.name.replaceAll(/\.(fasta|fa|fna)$/,'')
	  				[[sample_id:id],it]
	  		})
	  	}	  	
	  	if (opts.reads.illumina) {
				sr_ch = Channel.fromFilePairs(opts.reads.illumina,size:-1) { 
						file -> file.name.replaceAll(/_(R?[12])(_[0-9][0-9][0-9])?\.(fq|fastq)\.gz$/, '') 
				}
				.map({id,x -> [[sample_id:id],x]})
	  	}
			return fromChannels(sr_ch,lr_ch,fa_ch)			
		}
	}
	
}


