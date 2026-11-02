nextflow.enable.types = true

include { FLYE             } from './modules/flye'
include { MEDAKA_CONSENSUS } from './modules/medaka/consensus'
include { PILON_POLISH as PILON_POLISH_ROUND1 } from './subworkflows/pilon_polish'
include { PILON_POLISH as PILON_POLISH_ROUND2 } from './subworkflows/pilon_polish'
include { PILON_POLISH as PILON_POLISH_ROUND3 } from './subworkflows/pilon_polish'


// Short reads sample: read2 is null for single-end reads
record ShortReads {
	meta: Map
	reads: Path
	read2: Path?
}

// Long reads sample: one FASTQ per sample
record LongReads {
	meta: Map
	reads: Path
}

// Assembly of a sample: assembler output folder and assembly FASTA
record Assembly {
	meta: Map
	dir: Path
	fasta: Path
}

// A file or folder to put in the assembly folder under the given name
record NamedFile {
	name: String
	path: Path
}


// Gather the outputs of the assembly steps of a sample into a folder, each one under its given name
process ASSEMBLY_FOLDER {
    executor 'local'
    input:
    record(meta: Map, files: List<NamedFile>)
    stage:
    stageAs(files.collect({ f -> f.path }), 'in/?/*')
    output:
    record(meta: meta, dir: file('assembly'))
    script:
    def moves = files.collect({ f -> "mv '${f.path}' 'assembly/${f.name}'" }).join('\n')
    """
    mkdir assembly
    ${moves}
    """
}


workflow LONG_FLYE {
	take:
	long_reads: Channel<LongReads>
	flye_args: String?

	main:
	def flye_ch = FLYE(long_reads.map({ r -> r + record(args: flye_args) }))

	emit:
	assemblies: Channel<Assembly> = flye_ch
}


workflow LONG_FLYE_MEDAKA {
	take:
	long_reads: Channel<LongReads>
	flye_args: String?
	medaka_args: String?

	main:
	def flye_ch = LONG_FLYE(long_reads, flye_args)
	def medaka_ch = MEDAKA_CONSENSUS(
		flye_ch
			.map({ r -> record(meta: r.meta, fasta: r.fasta) })
			.join(long_reads, by: 'meta')
			.map({ r -> r + record(args: medaka_args) })
	)
	def folder_ch = ASSEMBLY_FOLDER(
		flye_ch.map({ r -> record(meta: r.meta, flye: r.dir) })
			.join(medaka_ch.map({ r -> record(meta: r.meta, medaka: r.dir) }), by: 'meta', remainder: true)
			.map({ r -> record(
				meta: r.meta,
				files: [record(name: '01_flye', path: r.flye)] +
					(r.medaka ? [record(name: '02_medaka', path: r.medaka)] : [])
			) })
	)

	emit:
	assemblies: Channel<Assembly> = medaka_ch
		.map({ r -> record(meta: r.meta, fasta: r.fasta) })
		.join(folder_ch, by: 'meta')
}


workflow HYBRID_FLYE_MEDAKA_PILON {
	take:
	short_reads: Channel<ShortReads>
	long_reads: Channel<LongReads>
	flye_args: String?
	medaka_args: String?

	main:
	def medaka_ch = LONG_FLYE_MEDAKA(long_reads, flye_args, medaka_args)
	def pilon1_ch = PILON_POLISH_ROUND1(medaka_ch, short_reads)
	def pilon2_ch = PILON_POLISH_ROUND2(pilon1_ch, short_reads)
	def pilon3_ch = PILON_POLISH_ROUND3(pilon2_ch, short_reads)

	def folder_ch = ASSEMBLY_FOLDER(
		medaka_ch.map({ r -> record(meta: r.meta, dir0: r.dir) })
			.join(pilon1_ch.map({ r -> record(meta: r.meta, dir1: r.dir) }), by: 'meta', remainder: true)
			.join(pilon2_ch.map({ r -> record(meta: r.meta, dir2: r.dir) }), by: 'meta', remainder: true)
			.join(pilon3_ch.map({ r -> record(meta: r.meta, dir3: r.dir) }), by: 'meta', remainder: true)
			.map({ r -> record(
				meta: r.meta,
				files: [record(name: '01_flye_medaka', path: r.dir0)] +
					(r.dir1 ? [record(name: '02_pilon_round1', path: r.dir1)] : []) +
					(r.dir2 ? [record(name: '03_pilon_round2', path: r.dir2)] : []) +
					(r.dir3 ? [record(name: '04_pilon_round3', path: r.dir3)] : [])
			) })
	)

	// last Pilon round, or Flye+Medaka assembly when polishing failed
	def fasta_ch = medaka_ch.map({ r -> record(meta: r.meta, fasta0: r.fasta) })
		.join(pilon3_ch.map({ r -> record(meta: r.meta, fasta3: r.fasta) }), by: 'meta', remainder: true)
		.map({ r -> record(meta: r.meta, fasta: r.fasta3 ?: r.fasta0) })

	emit:
	assemblies: Channel<Assembly> = fasta_ch.join(folder_ch, by: 'meta')
}
