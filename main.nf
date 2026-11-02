

include { SAMTOOLS_FASTQ as CONVERT_LONG_BAM_TO_FASTQ } from './modules/samtools/fastq' 
include { READS_FILTER; ReadsFilterParams             } from './wf/reads_filter/main.nf'
include { READS_QC as FILTERED_READS_QC               } from './wf/reads_qc/main.nf'
include { ASSEMBLE; AssemblerParams                   } from './wf/assemble/main.nf'
include { ASSEMBLIES_QC; AssembliesQcParams           } from './wf/assemblies_qc/main.nf'
include { AMR_ANNOT; AmrParams                       } from './wf/amr_annot/main.nf'
include { IDENTITY as COPY_ASSEMBLY_FASTA             } from './modules/identity/main.nf'
include { IDENTITY as COPY_SHORT_READS                } from './modules/identity/main.nf'
include { IDENTITY as COPY_LONG_READS                 } from './modules/identity/main.nf'

// ------------------------------------------------------------------
// Main entry point when running the pipeline from command line
// ------------------------------------------------------------------
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

record ReadsParams {
		ont: String?
    illumina: String?
}

params {
	verbose: Boolean = false
	csv:     Path? = null
	reads:   ReadsParams = record(ont:null, illumina:null)
	reads_filter: ReadsFilterParams = record(skip: false, long_limit: 1_000_000_000, short_limit: 1_000_000_000, chopper_args: '--quality 15 --minlength 200')
	assembler: AssemblerParams = record(name: 'long_flye_medaka', args: record(flye: null, medaka: '--bacteria'))
	assemblies_qc: AssembliesQcParams = record(checkm2: record(skip: false, db: null, args: ''))
	amr: AmrParams = record(
		orgfinder:     record(skip: false, args: ''),
		speciator:     record(skip: true, args: ''),
		amrfinderplus: record(skip: false, args: ''),
		mobtyper:      record(skip: false, args: ''),
		plasmidfinder: record(skip: false, skip_longread: false, skip_shortread: false, args: ''),
		resfinder:     record(skip: false, skip_longread: false, skip_shortread: false, args: ''),
		cgemlst:       record(skip: false, args: ''),
		MLST:          record(skip: false, args: ''),
		prokka:        record(skip: true, args: '--kingdom Bacteria')
	)
}

// Convert [meta, reads] tuple channels into the reads records expected by typed subworkflows
def toShortReads(ch) {
	// a single-end sample comes as a single Path, a paired-end one as [R1, R2]
	return ch.map({ m, f ->
		def fs = f instanceof List ? f : [f]
		record(meta: m, reads: fs[0], read2: fs.size() > 1 ? fs[1] : null)
	})
}

def toLongReads(ch) {
	return ch.map({ m, f -> record(meta: m, reads: f) })
}

// Assemblies are identified by a single meta map holding sample_id and assembler_name in typed subworkflows,
// and by a [[sample_id], [assembler_name]] pair in the output block
def toPairMeta(meta) {
	return [meta.subMap('sample_id'), meta.subMap('assembler_name')]
}

workflow {
	main:
		// Validate parameters and print summary of supplied ones
		validateParameters()
		log.info(paramsSummaryLog(workflow))

		// Parse parameters and get a channel of samples objects
		def samples = Samples.fromParams(params,{sheet,schema -> samplesheetToList(sheet, schema)})
		
		// Extract reads channels from samples channel
		def sr_ch = samples
				.filter({s -> s.reads_short_1})
				.map({s -> [s.subMap('sample_id'),[s.reads_short_1,s.reads_short_2].findAll({f -> f})]})
				| COPY_SHORT_READS
		
		def lr_ch = samples
			.filter({s -> s.reads_long})
			.map({s -> [s.subMap('sample_id'),s.reads_long]})
			| COPY_LONG_READS

		// CONVERT long_reads given in BAM/CRAM format into FASTQ format
		lr_ch = lr_ch.branch {_meta,f -> 
			bam: f.name =~ /\.(bam|cram)$/
			fq: true
		}
		lr_ch = lr_ch.fq.mix(CONVERT_LONG_BAM_TO_FASTQ(lr_ch.bam))

		// Reads as records for the typed subworkflows
		def sr_rec = toShortReads(sr_ch)
		def lr_rec = toLongReads(lr_ch)

		// Filter reads
		if (!params.reads_filter.skip) {
			def rd_flt_ch = READS_FILTER(params.reads_filter, sr_rec, lr_rec)
			sr_rec = rd_flt_ch.short_filtered
			lr_rec = rd_flt_ch.long_filtered
		}
		def reads_multiqc_html_ch = channel.empty()
		if (!params.reads_qc.skip) {
			def qc_ch = FILTERED_READS_QC(sr_rec, lr_rec)
			reads_multiqc_html_ch = qc_ch.multiqc_html
		}

		// Back to [meta, reads] tuples for the subworkflows not migrated yet
		// (a single-end sample is given as a single Path, as before)
		sr_ch = sr_rec.map({ r -> [r.meta, r.read2 ? [r.reads, r.read2] : r.reads] })
		lr_ch = lr_rec.map({ r -> [r.meta, r.reads] })

		
		// Run denovo assembly process when 
		// Determine assemblies that are done and thus that have to be run
		def asm_ch = samples.branch {s ->
			done: s.assembly_fasta
			todo: true
		}
		def asm_cpy = asm_ch.done
			.map({s -> [s.subMap('sample_id'),s.assembly_fasta]}) 
			| COPY_ASSEMBLY_FASTA
		def assemblies = ASSEMBLE(
				params.assembler,
				toShortReads(sr_ch.join(asm_ch.todo.map({s -> [s.subMap('sample_id')]}))),
				toLongReads(lr_ch.join(asm_ch.todo.map({s -> [s.subMap('sample_id')]})))
		)
		def asm_fa_ch = channel.empty().mix(
				asm_cpy.map({m,x -> [[m,[assembler_name:'none']],x]}),
				assemblies.map({r -> [[r.meta,[assembler_name:params.assembler.name]],r.fasta]})
		)
		def asm_dir_channel = assemblies.map({r -> [[r.meta,[assembler_name:params.assembler.name]],r.dir]})

		// Assemblies as records, with a single meta map holding sample_id and assembler_name
		def asm_rec = asm_fa_ch.map({m,x -> record(meta: m[0] + m[1], fasta: x)})

		// Compute assemblies QC stats
		def asm_qc_ch = ASSEMBLIES_QC(
			params.assemblies_qc,
			asm_rec,
			toShortReads(sr_ch.combine(asm_fa_ch.map({m,_x -> m}),by:0).map({m,r,a -> [m + a,r]})),
			toLongReads(lr_ch.combine(asm_fa_ch.map({m,_x -> m}),by:0).map({m,r,a -> [m + a,r]}))
		)

		// Finally run AMR annotations
		def amr_annot_ch = AMR_ANNOT(params.amr, asm_rec, sr_rec, lr_rec)
		
		

	publish:
		checkm2_db = asm_qc_ch.checkm2_db
		
		reads_long_filtered = lr_ch
		reads_short_filtered = sr_ch
		filtered_reads_multiqc_html = reads_multiqc_html_ch.ifEmpty([])
		
		assemblies_fasta        = asm_fa_ch
		assemblies_dir          = asm_dir_channel
		assemblies_checkm2      = asm_qc_ch.checkm2.map({r -> [toPairMeta(r.meta), r.report]})
		assemblies_multiqc_txt  = asm_qc_ch.assembly_multiqc_txt.ifEmpty([])
		assemblies_multiqc_xlsx = asm_qc_ch.assembly_multiqc_xlsx.ifEmpty([])

		resfinder            = amr_annot_ch.resfinder.map({r -> [toPairMeta(r.meta), r.dir]})
		resfinder_long       = amr_annot_ch.resfinder_long.map({r -> [r.meta, r.dir]})
		resfinder_short      = amr_annot_ch.resfinder_short.map({r -> [r.meta, r.dir]})
		plasmidfinder        = amr_annot_ch.plasmidfinder.map({r -> [toPairMeta(r.meta), r.dir]})
		plasmidfinder_long   = amr_annot_ch.plasmidfinder_long.map({r -> [r.meta, r.dir]})
		plasmidfinder_short  = amr_annot_ch.plasmidfinder_short.map({r -> [r.meta, r.dir]})
		amrfinderplus        = amr_annot_ch.amrfinderplus.map({r -> [toPairMeta(r.meta), r.dir]})
		mobtyper             = amr_annot_ch.mobtyper.map({r -> [toPairMeta(r.meta), r.tsv]})
		MLST                 = amr_annot_ch.MLST.map({r -> [toPairMeta(r.meta), r.txt]})
		prokka               = amr_annot_ch.prokka.map({r -> [toPairMeta(r.meta), r.dir]})
		cgemlst              = amr_annot_ch.cgemlst.map({r -> [toPairMeta(r.meta), r.dir]})
		orgfinder            = amr_annot_ch.orgfinder.map({r -> [toPairMeta(r.meta), r.dir]})
		amr_multireport_html = amr_annot_ch.multireport_html.ifEmpty([])
		amr_multireport_xlsx = amr_annot_ch.multireport_tables.map({r -> r.xlsx}).ifEmpty([])
		amr_multireport_assemblies_tsv  = amr_annot_ch.multireport_tables.map({r -> r.assemblies_tsv}).ifEmpty([])
		amr_multireport_contigs_tsv     = amr_annot_ch.multireport_tables.map({r -> r.contigs_tsv}).ifEmpty([])
		amr_multireport_resistances_tsv = amr_annot_ch.multireport_tables.map({r -> r.resistances_tsv}).ifEmpty([])
		amr_multireport_plasmidfinder_hits_tsv = amr_annot_ch.multireport_tables.map({r -> r.plasmidfinder_hits_tsv}).ifEmpty([])
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
		enabled params.verbose
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




