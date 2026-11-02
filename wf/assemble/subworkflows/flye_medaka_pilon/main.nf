
include { FLYE             } from './modules/flye'
include { FLYE_ADAPT       } from './modules/flye/adapt'
include { MEDAKA_CONSENSUS } from './modules/medaka/consensus'
include { MEDAKA_ADAPT     } from './modules/medaka/adapt'
include { ORGANIZE_FILES   } from './modules/organize_files'
include { PILON_POLISH as PILON_POLISH_ROUND1 } from './subworkflows/pilon_polish'
include { PILON_POLISH as PILON_POLISH_ROUND2 } from './subworkflows/pilon_polish'
include { PILON_POLISH as PILON_POLISH_ROUND3 } from './subworkflows/pilon_polish'


workflow LONG_FLYE {
  take:
		fql_ch
	main:
		fql_ch
		| FLYE
		| FLYE_ADAPT
	emit:
		fasta = FLYE_ADAPT.out.fasta
		dir = FLYE.out
}

workflow LONG_FLYE_MEDAKA {
  take:
		fql_ch
	main:
		LONG_FLYE(fql_ch)
		LONG_FLYE.out.fasta.join(fql_ch)
		| MEDAKA_CONSENSUS 
		| MEDAKA_ADAPT
		
		ORGANIZE_FILES(
			LONG_FLYE.out.dir
			.join(MEDAKA_CONSENSUS.out, remainder: true)
			.map({meta,x1,x2 -> tuple(meta,[[x1,'01_flye'],[x2,'02_medaka']].findAll({x,_y -> x}))})
		)
	emit:
		fasta = MEDAKA_ADAPT.out.fasta
		dir = ORGANIZE_FILES.out		
}


workflow HYBRID_FLYE_MEDAKA_PILON {
	take:
		fqs_ch
		fql_ch
	main:
		LONG_FLYE_MEDAKA(fql_ch)
		PILON_POLISH_ROUND1(LONG_FLYE_MEDAKA.out.fasta,fqs_ch)
		PILON_POLISH_ROUND2(PILON_POLISH_ROUND1.out.fasta,fqs_ch)
		PILON_POLISH_ROUND3(PILON_POLISH_ROUND2.out.fasta,fqs_ch)
		ORGANIZE_FILES(
			LONG_FLYE_MEDAKA.out.dir
			.join(PILON_POLISH_ROUND1.out.dir, remainder: true)
			.join(PILON_POLISH_ROUND2.out.dir, remainder: true)
			.join(PILON_POLISH_ROUND3.out.dir, remainder: true)
			.map({meta,x1,x2,x3,x4 -> tuple(meta,[[x1,'01_flye_medaka'],[x2,'02_pilon_round1'],[x3,'03_pilon_round2'],[x4,'04_pilon_round3']].findAll({x,_y -> x}))})
		)
		fa_ch = LONG_FLYE_MEDAKA.out.fasta
		.join(PILON_POLISH_ROUND3.out.fasta,remainder: true)
		.map({meta,fa1,fa2 -> [meta,fa2?fa2:fa1]})
	emit:
		fasta = fa_ch
		dir = ORGANIZE_FILES.out
}
