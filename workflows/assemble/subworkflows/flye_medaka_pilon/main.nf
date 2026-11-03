nextflow.enable.types = true

include { FLYE             } from './modules/flye'
include { MEDAKA_CONSENSUS } from './modules/medaka/consensus'
include { PILON_POLISH as PILON_POLISH_ROUND1 } from './subworkflows/pilon_polish'
include { PILON_POLISH as PILON_POLISH_ROUND2 } from './subworkflows/pilon_polish'
include { PILON_POLISH as PILON_POLISH_ROUND3 } from './subworkflows/pilon_polish'


// Illumina reads of a sample: read2 is null for single-end reads
record ShortReads {
	meta: Map?
	read1: Path
	read2: Path?
}

// Arguments of the assembly tools (each assembler reads its own field)
record AssemblerArgs {
	flye: String?
	medaka: String?
	spades: String?
	unicycler: String?
	hybracter: String?
}

// Input of an assembler: the reads of a sample (short_reads2 is null for single-end reads) and the tools arguments
record AssemblyInput {
	meta: Map?
	short_reads: Path?
	short_reads2: Path?
	long_reads: Path?
	args: AssemblerArgs?
}

// Assembly of a sample: assembler output folder and assembly FASTA
record Genome {
	meta: Map
	dir: Path
	fasta: Path
	gff: Path?
	info: Path?
}

// A file or folder to put in the assembly folder under the given name
record NamedFile {
	name: String
	path: Path
}


// Gather the outputs of the assembly steps of a sample into a folder, each one under its given name
process ASSEMBLY_FOLDER {
    label 'local'
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
	samples: Channel<AssemblyInput>    // with long_reads

	main:
	def flye_ch = FLYE(samples.map({ r -> record(meta: r.meta, reads: r.long_reads, args: r.args?.flye) }))

	emit:
	assemblies: Channel<Genome> = flye_ch
}


workflow LONG_FLYE_MEDAKA {
	take:
	samples: Channel<AssemblyInput>    // with long_reads

	main:
	def flye_ch = LONG_FLYE(samples)
	def medaka_ch = MEDAKA_CONSENSUS(
		flye_ch
			.map({ r -> record(meta: r.meta, fasta: r.fasta) })
			.join(samples.map({ r -> record(meta: r.meta, reads: r.long_reads, args: r.args?.medaka) }), by: 'meta')
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
	assemblies: Channel<Genome> = medaka_ch
		.map({ r -> record(meta: r.meta, fasta: r.fasta) })
		.join(folder_ch, by: 'meta')
}


workflow HYBRID_FLYE_MEDAKA_PILON {
	take:
	samples: Channel<AssemblyInput>    // with long_reads, and short_reads for the Pilon polishing

	main:
	def short_reads = samples
		.filter({ r -> r.short_reads })
		.map({ r -> record(meta: r.meta, read1: r.short_reads, read2: r.short_reads2) })
	def medaka_ch = LONG_FLYE_MEDAKA(samples)
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
	assemblies: Channel<Genome> = fasta_ch.join(folder_ch, by: 'meta')
}
