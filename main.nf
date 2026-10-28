
include { SAMTOOLS_FASTQ as CONVERT_LONG_BAM_TO_FASTQ } from './modules/samtools/fastq' 
include { READS_FILTER                                } from './wf/reads_filter/main.nf'
include { READS_QC as FILTERED_READS_QC               } from './wf/reads_qc/main.nf'
include { ASSEMBLE                                    } from './wf/assemble/main.nf'
include { ASSEMBLIES_QC                               } from './wf/assemblies_qc/main.nf'
include { AMR_ANNOT                                   } from './wf/amr_annot/main.nf'
include { IDENTITY as COPY_ASSEMBLY_FASTA             } from './modules/identity/main.nf'
include { IDENTITY as COPY_SHORT_READS                } from './modules/identity/main.nf'
include { IDENTITY as COPY_LONG_READS                 } from './modules/identity/main.nf'

// ------------------------------------------------------------------
// Main entry point when running the pipeline from command line
// ------------------------------------------------------------------
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'


workflow {
	main:
		// Validate parameters and print summary of supplied ones
		validateParameters()
		log.info(paramsSummaryLog(workflow))

		// Parse parameters and get a channel of samples objects
		def samples = Samples.fromParams(params,{sheet,schema -> samplesheetToList(sheet, schema)})
		
		// Extract reads channels from samples channel
		def sr_ch = samples
				.filter({it.reads_short_1})
				.map({[it.subMap('sample_id'),[it.reads_short_1,it.reads_short_2].findAll({it})]})
				| COPY_SHORT_READS
		
		def lr_ch = samples
			.filter({it.reads_long})
			.map({[it.subMap('sample_id'),it.reads_long]})
			| COPY_LONG_READS

		// CONVERT long_reads given in BAM/CRAM format into FASTQ format
		lr_ch = lr_ch.branch {meta,f -> 
			bam: f.name =~ /\.(bam|cram)$/
			fq: true
		}
		lr_ch = lr_ch.fq.mix(CONVERT_LONG_BAM_TO_FASTQ(lr_ch.bam))

		// Filter reads
		if (!params.reads_filter.skip) {
			READS_FILTER(sr_ch,lr_ch)
			sr_ch = READS_FILTER.out.short_filtered
			lr_ch = READS_FILTER.out.long_filtered
		}
		def reads_multiqc_html_ch = Channel.empty()
		if (!params.reads_qc.skip) {
			FILTERED_READS_QC(sr_ch,lr_ch)
			reads_multiqc_html_ch = FILTERED_READS_QC.out.multiqc_html
		}
		
		
		// Run denovo assembly process when 
		// Determine assemblies that are done and thus that have to be run
		def asm_ch = samples.branch {
			done: it.assembly_fasta
			todo: true
		}
		asm_ch.done.map({[it.subMap('sample_id'),it.assembly_fasta]}) 
		| COPY_ASSEMBLY_FASTA
		ASSEMBLE(
				params.assembler.name,
				sr_ch.join(asm_ch.todo.map({[it.subMap('sample_id')]})),
				lr_ch.join(asm_ch.todo.map({[it.subMap('sample_id')]}))
		)
		def asm_fa_ch = Channel.empty().mix(
				COPY_ASSEMBLY_FASTA.out.map({[[it[0],[assembler_name:'none']],it[1]]}),
				ASSEMBLE.out.fasta.map({[[it[0],[assembler_name:params.assembler.name]],it[1]]})
		)
		def asm_dir_channel = ASSEMBLE.out.dir.map({[[it[0],[assembler_name:params.assembler.name]],it[1]]})

		// Compute assemblies QC stats
		ASSEMBLIES_QC(
			params.assemblies_qc,
			asm_fa_ch,
			sr_ch.combine(asm_fa_ch.map({m,x->m}),by:0).map({[[it[0],it[2]],it[1]]}),
			lr_ch.combine(asm_fa_ch.map({m,x->m}),by:0).map({[[it[0],it[2]],it[1]]})
		)

		// Finally run AMR annotations
		AMR_ANNOT(params.amr,asm_fa_ch,sr_ch,lr_ch)
		
		

	publish:
		checkm2_db = ASSEMBLIES_QC.out.checkm2_db
		
		reads_long_filtered = lr_ch
		reads_short_filtered = sr_ch
		filtered_reads_multiqc_html = reads_multiqc_html_ch.ifEmpty([])
		
		assemblies_fasta        = asm_fa_ch
		assemblies_dir          = asm_dir_channel
		assemblies_checkm2      = ASSEMBLIES_QC.out.checkm2
		assemblies_multiqc_txt  = ASSEMBLIES_QC.out.assembly_multiqc_txt.ifEmpty([])
		assemblies_multiqc_xlsx = ASSEMBLIES_QC.out.assembly_multiqc_xlsx.ifEmpty([])

		resfinder            = AMR_ANNOT.out.resfinder
		resfinder_long       = AMR_ANNOT.out.resfinder_long
		resfinder_short      = AMR_ANNOT.out.resfinder_short
		plasmidfinder        = AMR_ANNOT.out.plasmidfinder
		plasmidfinder_long   = AMR_ANNOT.out.plasmidfinder_long
		plasmidfinder_short  = AMR_ANNOT.out.plasmidfinder_short
		amrfinderplus        = AMR_ANNOT.out.amrfinderplus
		mobtyper             = AMR_ANNOT.out.mobtyper
		MLST                 = AMR_ANNOT.out.MLST
		prokka               = AMR_ANNOT.out.prokka
		cgemlst              = AMR_ANNOT.out.cgemlst
		orgfinder            = AMR_ANNOT.out.orgfinder
		amr_multireport_html = AMR_ANNOT.out.multireport_html.ifEmpty([])
		amr_multireport_xlsx = AMR_ANNOT.out.multireport_xlsx.ifEmpty([])
		amr_multireport_assemblies_tsv  = AMR_ANNOT.out.multireport_assemblies_tsv.ifEmpty([])
		amr_multireport_contigs_tsv     = AMR_ANNOT.out.multireport_contigs_tsv.ifEmpty([])
		amr_multireport_resistances_tsv = AMR_ANNOT.out.multireport_resistances_tsv.ifEmpty([])
		amr_multireport_plasmidfinder_hits_tsv = AMR_ANNOT.out.multireport_plasmidfinder_hits_tsv.ifEmpty([])
}

output {
	checkm2_db {
		path { x -> x >> "db/checkm2"}
	}
	
	reads_long_filtered {
		path { m,x -> x >> "reads/filtered_long/${m.sample_id}.fastq.gz"}
	}
	reads_short_filtered {
		path { m,x -> 
			x = x instanceof List?x:[x]
			x.eachWithIndex { xi, i ->
      	xi >> "reads/filtered_short/${m.sample_id}_R${i+1}.fastq.gz"
      }
		}
	}

	filtered_reads_multiqc_html {
		path { "reports/" }
	}

	assemblies_fasta {
		mode "copy" //because link is problematic
		path { m,x -> x >> "assemblies/${m[0].sample_id}__${m[1].assembler_name}.fasta"}
	}
	assemblies_dir {
		mode "copy" //because link is problematic
		enabled params.verbosity>0
		path { m,x -> x >> "assemblies/${m[0].sample_id}__${m[1].assembler_name}.output"}
	}	
	assemblies_checkm2 {
		path { m,x -> x >> "assemblies/${m[0].sample_id}__${m[1].assembler_name}.qc/checkm2_quality.tsv"}
	}
	assemblies_multiqc_txt {
		path { "reports/"}
	}
	assemblies_multiqc_xlsx {
		path { "reports/"}
	}

	
	resfinder_long {
		path { m,x -> x >> "amr/reads/${m.sample_id}/resfinder_long"}
	}
	resfinder_short {
		path { m,x -> x >> "amr/reads/${m.sample_id}/resfinder_short"}
	}
	plasmidfinder_long {
		path { m,x -> x >> "amr/reads/${m.sample_id}/plasmidfinder_long"}
	}
	plasmidfinder_short {
		path { m,x -> x >> "amr/reads/${m.sample_id}/plasmidfinder_short"}
	}
	

	resfinder {
		path { m,x -> x >> "amr/assemblies/${m[0].sample_id}__${m[1].assembler_name}.amr/resfinder"}
	}
	plasmidfinder {
		path { m,x -> x >> "amr/assemblies/${m[0].sample_id}__${m[1].assembler_name}.amr/plasmidfinder"}
	}
	amrfinderplus {
		path { m,x -> x >> "amr/assemblies/${m[0].sample_id}__${m[1].assembler_name}.amr/amrfinderplus"}
	}
	mobtyper {
		path { m,x -> x >> "amr/assemblies/${m[0].sample_id}__${m[1].assembler_name}.amr/mobtyper"}
	}
	MLST {
		path { m,x -> x >> "amr/assemblies/${m[0].sample_id}__${m[1].assembler_name}.amr/MLST"}
	}
	prokka {
		path { m,x -> x >> "amr/assemblies/${m[0].sample_id}__${m[1].assembler_name}.amr/prokka"}
	}
	orgfinder {
		path { m,x -> x >> "amr/assemblies/${m[0].sample_id}__${m[1].assembler_name}.amr/orgfinder"}
	}
	cgemlst {
		path { m,x -> x >> "amr/assemblies/${m[0].sample_id}__${m[1].assembler_name}.amr/cgemlst"}
	}
	
	amr_multireport_html {
		path { "reports/"}
	}
	amr_multireport_xlsx {
		path { "reports/"}
	}
	amr_multireport_assemblies_tsv {
		path { "reports/"}
	}	
	amr_multireport_contigs_tsv {
		path { "reports/"}
	}	
	amr_multireport_resistances_tsv {
		path { "reports/"}
	}	
  amr_multireport_plasmidfinder_hits_tsv {
		path { "reports/"}
	}
}




