

include { LONG_FLYE                     } from './subworkflows/flye_medaka_pilon'
include { LONG_FLYE_MEDAKA              } from './subworkflows/flye_medaka_pilon'
include { UNICYCLER as LONG_UNICYCLER   } from './subworkflows/unicycler'
include { HYBRACTER as LONG_HYBRACTER   } from './subworkflows/hybracter'
include { UNICYCLER as SHORT_UNICYCLER  } from './subworkflows/unicycler'
include { SPADES    as SHORT_SPADES     } from './subworkflows/spades'
include { UNICYCLER as HYBRID_UNICYCLER } from './subworkflows/unicycler'
include { HYBRACTER as HYBRID_HYBRACTER } from './subworkflows/hybracter'
include { HYBRID_FLYE_MEDAKA_PILON      } from './subworkflows/flye_medaka_pilon'


workflow ASSEMBLE {
	take:
		assembler_name
		fqs_ch    // channel: [ val(meta), path(short_reads) ]
		fql_ch    // channel: [ val(meta), path(long_reads) ]
	main:
		def assemblies_fasta = channel.empty()
		def assemblies_dir   = channel.empty()
		
		if (assembler_name == 'long_flye') {
			LONG_FLYE(fql_ch)
			assemblies_fasta = LONG_FLYE.out.fasta
			assemblies_dir   = LONG_FLYE.out.dir
		} else if (assembler_name == 'long_flye_medaka') {
			LONG_FLYE_MEDAKA(fql_ch)
			assemblies_fasta = LONG_FLYE_MEDAKA.out.fasta
			assemblies_dir   = LONG_FLYE_MEDAKA.out.dir
		} else if (assembler_name == 'long_unicycler') {
			LONG_UNICYCLER(channel.empty(),fql_ch)
			assemblies_fasta = LONG_UNICYCLER.out.fasta
			assemblies_dir   = LONG_UNICYCLER.out.dir
		} else if (assembler_name == 'long_hybracter') {
			LONG_HYBRACTER(channel.empty(),fql_ch)
			assemblies_fasta = LONG_HYBRACTER.out.fasta
			assemblies_dir   = LONG_HYBRACTER.out.dir
		} else if (assembler_name == 'short_unicycler') {
			SHORT_UNICYCLER(fqs_ch,channel.empty())
			assemblies_fasta = SHORT_UNICYCLER.out.fasta
			assemblies_dir   = SHORT_UNICYCLER.out.dir
		} else if (assembler_name == 'short_spades') {
			SHORT_SPADES(fqs_ch,channel.empty())
			assemblies_fasta = SHORT_SPADES.out.fasta
			assemblies_dir   = SHORT_SPADES.out.dir
		} else if (assembler_name == 'hybrid_unicycler') {
			HYBRID_UNICYCLER(fqs_ch,fql_ch)
			assemblies_fasta = HYBRID_UNICYCLER.out.fasta
			assemblies_dir   = HYBRID_UNICYCLER.out.dir
		} else if (assembler_name == 'hybrid_hybracter') {
			HYBRID_HYBRACTER(fqs_ch,fql_ch)
			assemblies_fasta = HYBRID_HYBRACTER.out.fasta
			assemblies_dir   = HYBRID_HYBRACTER.out.dir
		} else if (assembler_name == 'hybrid_flye_medaka_pilon') {
			HYBRID_FLYE_MEDAKA_PILON(fqs_ch,fql_ch)
			assemblies_fasta = HYBRID_FLYE_MEDAKA_PILON.out.fasta
			assemblies_dir   = HYBRID_FLYE_MEDAKA_PILON.out.dir
		} else {
			error "Unknown assembler name :${assembler_name}"
		}
		
	emit:
		fasta = assemblies_fasta
		dir = assemblies_dir
}

