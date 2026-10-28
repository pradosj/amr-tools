#!/usr/bin/env nextflow

include { GZIP_DECOMPRESS as GZIP_DECOMPRESS_FASTA } from './modules/gzip'
include { RESFINDER                                } from './modules/cgetools/resfinder'
include { RESFINDER     as RESFINDER_LONGREAD      } from './modules/cgetools/resfinder'
include { RESFINDER     as RESFINDER_SHORTREAD     } from './modules/cgetools/resfinder'
include { PLASMIDFINDER                            } from './modules/cgetools/plasmidfinder'
include { PLASMIDFINDER as PLASMIDFINDER_LONGREAD  } from './modules/cgetools/plasmidfinder'
include { PLASMIDFINDER as PLASMIDFINDER_SHORTREAD } from './modules/cgetools/plasmidfinder'
include { ORGFINDER_DETECT                         } from './modules/orgfinder/detect'
include { SPECIATOR                                } from './modules/speciator'
include { AMRFINDERPLUS_UPDATE                     } from './modules/amrfinderplus/update'
include { AMRFINDERPLUS_RUN                        } from './modules/amrfinderplus/run'
include { PROKKA                                   } from './modules/tseemann/prokka'
include { MLST                                     } from './modules/tseemann/mlst'
include { MLST as CGEMLST                          } from './modules/cgetools/mlst'
include { MOBTYPER                                 } from './modules/mobsuite/mobtyper.nf'
include { SAMTOOLS_FAIDX                           } from './modules/samtools/faidx'
include { MULTIREPORT                              } from './subworkflows/multireport'




workflow AMR_ANNOT_ASSEMBLY {
	take:
		opts
		asm_ch
	main:

		def orgArgs = new groovy.json.JsonSlurper().parseText(file("${moduleDir}/assets/default_org_args.json").text)

		// CGE - RESFINDER
		def resfinder_ch = RESFINDER(asm_ch.filter({!opts.resfinder.skip}),'fasta')

		// Plasmid typing
		def plasmidfinder_ch = asm_ch.filter({!opts.plasmidfinder.skip}) | PLASMIDFINDER

		// NCBI AMRfinder+
		def amrfinderplus_ch = Channel.empty()
		if (!opts.amrfinderplus.skip) {
			def amrfinderplus_db = AMRFINDERPLUS_UPDATE()
			amrfinderplus_ch = AMRFINDERPLUS_RUN(asm_ch,amrfinderplus_db)
		}

		// MOBsuite - MOBtyper
		def mobtyper_ch = Channel.empty()
		if (!opts.mobtyper.skip) {
			mobtyper_ch = asm_ch | MOBTYPER
		}

		// PROKKA annotations
		def prokka_ch = Channel.empty()
		if (!opts.prokka.skip) {
			prokka_ch = asm_ch | PROKKA
		}

		// Speciator
		def speciator_ch = Channel.empty()
		if (!opts.speciator.skip) {
			speciator_ch = asm_ch | SPECIATOR	
		}
		
		def MLST_ch = Channel.empty()
		if (!opts.MLST.skip) {
			MLST_ch = asm_ch | MLST
		}

    // Run orgfinder to auto detect organism
    def orgfinder_dir_ch = Channel.empty()
    def orgfinder_name_ch = Channel.empty()
    if (!opts.orgfinder.skip) {
    	ORGFINDER_DETECT(asm_ch)
    	orgfinder_dir_ch = ORGFINDER_DETECT.out.orgfinder
    	orgfinder_name_ch = ORGFINDER_DETECT.out.org_name
    }

		// MLST typing on detected org_name
		def cgemlst_ch = Channel.empty()
		if (!opts.cgemlst.skip) {
			cgemlst_ch = asm_ch
				.join(orgfinder_name_ch)
				.map({m,fa,org -> def arg=orgArgs[org]?:[:];[m,fa,arg.cgemlst]})
				.filter({m,fa,arg -> arg})
				| CGEMLST
		}

	emit:
			orgfinder           = orgfinder_dir_ch
			amrfinderplus       = amrfinderplus_ch
			resfinder           = resfinder_ch
			mobtyper            = mobtyper_ch
			plasmidfinder       = plasmidfinder_ch
			prokka              = prokka_ch
			org_name            = orgfinder_name_ch
			cgemlst             = cgemlst_ch
			MLST                = MLST_ch
}



workflow AMR_ANNOT_READS {
	take:
		opts
		fqs_ch
		fql_ch
	main:
		PLASMIDFINDER_LONGREAD(fql_ch.filter({!opts.plasmidfinder.skip_longread}))
		RESFINDER_LONGREAD(fql_ch.filter({!opts.resfinder.skip_longread}),'nanopore')
		PLASMIDFINDER_SHORTREAD(fqs_ch.filter({!opts.plasmidfinder.skip_shortread}))
		RESFINDER_SHORTREAD(fqs_ch.filter({!opts.resfinder.skip_shortread}),'illumina')
	emit:
			resfinder_long      = RESFINDER_LONGREAD.out
			resfinder_short     = RESFINDER_SHORTREAD.out
			plasmidfinder_long  = PLASMIDFINDER_LONGREAD.out
			plasmidfinder_short = PLASMIDFINDER_SHORTREAD.out
}




workflow AMR_ANNOT {
	take:
		opts
		asm_ch   // tuple(meta,path()) with meta=[[sample_id:""],[assembler_name:""]]
		fqs_ch   // tuple(meta,path()) with meta=[sample_id:""]
		fql_ch   // tuple(meta,path()) with meta=[sample_id:""]
	main:
		fai_ch = SAMTOOLS_FAIDX(asm_ch)
		AMR_ANNOT_ASSEMBLY(opts,asm_ch)
		AMR_ANNOT_READS(opts,fqs_ch,fql_ch)
		MULTIREPORT(
			asm_ch,
			fai_ch,
			AMR_ANNOT_ASSEMBLY.out.org_name,
			AMR_ANNOT_ASSEMBLY.out.orgfinder,
			AMR_ANNOT_ASSEMBLY.out.amrfinderplus,
			AMR_ANNOT_ASSEMBLY.out.resfinder,
			AMR_ANNOT_ASSEMBLY.out.mobtyper,
			AMR_ANNOT_ASSEMBLY.out.plasmidfinder,
			AMR_ANNOT_ASSEMBLY.out.cgemlst,
			AMR_ANNOT_ASSEMBLY.out.MLST,
			AMR_ANNOT_ASSEMBLY.out.prokka,
			AMR_ANNOT_READS.out.resfinder_long,
			AMR_ANNOT_READS.out.plasmidfinder_long,
			AMR_ANNOT_READS.out.resfinder_short,
			AMR_ANNOT_READS.out.plasmidfinder_short
		)
	emit:
			orgfinder           = AMR_ANNOT_ASSEMBLY.out.orgfinder
			amrfinderplus       = AMR_ANNOT_ASSEMBLY.out.amrfinderplus
			resfinder           = AMR_ANNOT_ASSEMBLY.out.resfinder
			mobtyper            = AMR_ANNOT_ASSEMBLY.out.mobtyper
			plasmidfinder       = AMR_ANNOT_ASSEMBLY.out.plasmidfinder
			cgemlst             = AMR_ANNOT_ASSEMBLY.out.cgemlst
			MLST                = AMR_ANNOT_ASSEMBLY.out.MLST
			prokka              = AMR_ANNOT_ASSEMBLY.out.prokka
			resfinder_long      = AMR_ANNOT_READS.out.resfinder_long
			resfinder_short     = AMR_ANNOT_READS.out.resfinder_short
			plasmidfinder_long  = AMR_ANNOT_READS.out.plasmidfinder_long
			plasmidfinder_short = AMR_ANNOT_READS.out.plasmidfinder_short

	    multireport_folder = MULTIREPORT.out.folder
    	multireport_html   = MULTIREPORT.out.html.map({m,x -> x})
    	multireport_xlsx   = MULTIREPORT.out.xlsx.map({m,x -> x})
    	multireport_contigs_tsv = MULTIREPORT.out.contigs_tsv.map({m,x -> x})
    	multireport_resistances_tsv = MULTIREPORT.out.resistances_tsv.map({m,x -> x})
    	multireport_assemblies_tsv = MULTIREPORT.out.assemblies_tsv.map({m,x -> x})
}














