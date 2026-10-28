
include { RMD_RENDER        } from '../../modules/rmd/render'
include { ORGANIZE_FILES    } from '../../modules/organize_files'

process MULTITABLE {
	container 'registry.gitlab.unige.ch/amr-genomics/rscript:main'
    memory '8 GB'
    cpus 2
    time '30 min'
    input:
		tuple val(meta),path("db")
		each path("lib_typing.R")
    output:
        tuple val(meta), path('multitable.xlsx'), emit:'xlsx'
        tuple val(meta), path('multitable_contigs.tsv'), emit:'contigs_tsv'
        tuple val(meta), path('multitable_resistances.tsv'), emit:'resistances_tsv'
    script:
        """
		#!/usr/bin/env Rscript
		source("lib_typing.R")
		db <- db_load("db")
		tbl <- list(
			assemblies = summarise_assembly(db),
			contigs = summarise_contigs(db),
			resistances = summarise_resistances(db),
			plasmidfinder_hits = summarise_plasmidfinder_hits(db)
		)
		tbl\$contigs <- left_join(tbl\$contigs,select(tbl\$assemblies,sample_id,assembler_name,species_name,mlst_type),by=c('sample_id','assembler_name'))
		openxlsx::write.xlsx(tbl,file="multitable.xlsx")
		readr::write_tsv(tbl\$contigs,file="multitable_contigs.tsv")
		readr::write_tsv(tbl\$contigs,file="multitable_resistances.tsv")
        """
}


workflow MULTIREPORT {
	take:
		fa
		fai
		org_name
		orgfinder
		amrfinderplus
		resfinder
		mobtyper
		plasmidfinder
		cgemlst
		MLST
		prokka
		long_resfinder
		long_plasmidfinder
		short_resfinder
		short_plasmidfinder
		
	main:
		ORGANIZE_FILES(
			Channel.empty().mix(
				fa.map(           {m,file -> [file,"amr/assemblies/${m[0].sample_id}__${m[1].assembler_name}.fasta"]}),
				fai.map(          {m,file -> [file,"amr/assemblies/${m[0].sample_id}__${m[1].assembler_name}.fasta.fai"]}),
				org_name.map(     {m,val  -> [val ,"amr/assemblies/${m[0].sample_id}__${m[1].assembler_name}.amr/org_name.txt"]}),
				orgfinder.map(    {m,file -> [file,"amr/assemblies/${m[0].sample_id}__${m[1].assembler_name}.amr/${file.name}"]}),
				amrfinderplus.map({m,file -> [file,"amr/assemblies/${m[0].sample_id}__${m[1].assembler_name}.amr/${file.name}"]}),
				resfinder.map(    {m,file -> [file,"amr/assemblies/${m[0].sample_id}__${m[1].assembler_name}.amr/${file.name}"]}),
				plasmidfinder.map({m,file -> [file,"amr/assemblies/${m[0].sample_id}__${m[1].assembler_name}.amr/${file.name}"]}),
				cgemlst.map(      {m,file -> [file,"amr/assemblies/${m[0].sample_id}__${m[1].assembler_name}.amr/${file.name}"]}),
				mobtyper.map(     {m,file -> [file,"amr/assemblies/${m[0].sample_id}__${m[1].assembler_name}.amr/${file.name}"]}),
				MLST.map(         {m,file -> [file,"amr/assemblies/${m[0].sample_id}__${m[1].assembler_name}.amr/${file.name}"]}),
				prokka.map(       {m,file -> [file,"amr/assemblies/${m[0].sample_id}__${m[1].assembler_name}.amr/${file.name}"]}),
				long_resfinder.map(     {m,file -> [file,"amr/reads/${m.sample_id}.reads_amr/resfinder_long"]}),
				long_plasmidfinder.map( {m,file -> [file,"amr/reads/${m.sample_id}.reads_amr/plasmidfinder_long"]}),
				short_resfinder.map(    {m,file -> [file,"amr/reads/${m.sample_id}.reads_amr/resfinder_short"]}),
				short_plasmidfinder.map({m,file -> [file,"amr/reads/${m.sample_id}.reads_amr/plasmidfinder_short"]})
			)
			.collect({[it]})
			.map({["aggregated_folder",it]})
		)

		RMD_RENDER(
			ORGANIZE_FILES.out.map({m,x -> ["multireport.html",x,"indir='${x}'"]}),
			file("${moduleDir}/assets/multireport.Rmd"),
			file("${moduleDir}/assets/lib_typing.R")
		)

		MULTITABLE(
			ORGANIZE_FILES.out.map({m,x -> ["multitable.xlsx",x]}),
			file("${moduleDir}/assets/lib_typing.R")
		)

	emit:
		folder = ORGANIZE_FILES.out
		html   = RMD_RENDER.out       // channel: [ meta, path(html) ]
		xlsx   = MULTITABLE.out.xlsx  // channel: [ meta, path(xlsx) ]
		contigs_tsv     = MULTITABLE.out.contigs_tsv      // channel: [ meta, path(tsv) ]
		resistances_tsv = MULTITABLE.out.resistances_tsv  // channel: [ meta, path(tsv) ]
}
