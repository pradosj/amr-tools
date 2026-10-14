
include { NCBI_DATASET_DOWNLOAD_GENOME } from './modules/ncbi/dataset/main.nf'
include { NCBI_TAXDUMP_DOWNLOAD        } from './modules/ncbi/taxdump/main.nf'

// Make the tsv file with all accession numbers
process STEP1_ACCESSIONS_TSV {
    container 'docker.io/staphb/ncbi-datasets:18.18.0'
    memory '8 GB'
    cpus 1
    time '30 min'
    input:
  		path('genomes/query*')
  	output:
  		tuple(path('genomes/',type: 'dir'), path('accessions.tsv'))
    script:
	    """
			cat genomes/query*/ncbi_dataset/data/assembly_data_report.jsonl \
			  | dataformat tsv genome --force --fields accession,organism-name,organism-tax-id,assmstats-total-sequence-len,assmstats-total-number-of-chromosomes \
			  | sort \
			  | uniq \
			  > accessions.tsv
	    """
}

process STEP2_MOLECULE_TYPES_TSV {
		label 'rscript'
		container 'registry.gitlab.unige.ch/amr-genomics/rscript:v3'
    memory '2 GB'
    cpus 1
    time '30 min'
    input:
  		tuple(path('ncbi_db/genomes'),path('ncbi_db/accessions.tsv'))
  	output:
  		path('ncbi_db/',type: 'dir')
    script:
	    """
			jq -r '[.assemblyAccession, .chrName, .assignedMoleculeLocationType] | @tsv' ncbi_db/genomes/query*/ncbi_dataset/data/*/sequence_report.jsonl \
			| sort \
			| uniq \
			> ncbi_db/molecule_types.tsv
	    """
}

process STEP3_DB_BUILD {
	  label 'rscript'
	  container 'registry.gitlab.unige.ch/amr-genomics/rscript:v3'
    memory '8 GB'
    cpus 1
    time '30 min'
    input:
        path('ncbi_db')
        each path('taxdump')
        each path('db_build.R')
    output:
        path('orgfinder_db',type: 'dir', emit: 'dir')
        path('orgfinder_db/fna/*.fna', emit: 'fasta')
        path('orgfinder_db/db.tsv', emit: 'tsv')
    script:
				"""
				Rscript db_build.R
				"""
}


workflow ORGFINDER_DB {
	main:
		def taxdump = NCBI_TAXDUMP_DOWNLOAD()
		def genomes_ch = Channel.of(
			"taxon 'Pseudomonas aeruginosa'  --reference --assembly-level complete --include genome,seq-report",
			"taxon 'Acinetobacter baumannii' --reference --assembly-level complete --include genome,seq-report",
			"taxon 'Enterococcus'            --reference --assembly-level complete --include genome,seq-report",
			"taxon 'Staphylococcus'          --reference --assembly-level complete --include genome,seq-report",
			"taxon 'Streptococcus'           --reference --assembly-level complete --include genome,seq-report",
			"taxon 'Enterobacterales'        --reference --assembly-level complete --include genome,seq-report",
			"taxon 'Aeromonas'               --reference --assembly-level complete --include genome,seq-report",
			"taxon 'Myroides'                --reference --assembly-level complete --include genome,seq-report",
			"taxon 'Enterococcus faecalis'   --reference --include genome,seq-report",
			"taxon 'Citrobacter murliniae'   --reference --include genome,seq-report",
			'accession GCA_040096145.1', //E. intestinihominis	3133180	GCA_040096145.1	Genome of type strain CLA-AC-H004ᵀ
			'accession GCA_001875655.1', //E. hormaechei subsp. hormaechei	301105 / parent 158836	GCA_001875655.1	Genome of species type strain ATCC 49162ᵀ
			'accession GCA_001729745.1', //E. hormaechei subsp. hoffmannii	1812934	GCA_001729745.1	Complete genome of type strain DSM 14563ᵀ
			'accession GCF_048568405.1', //E. intestinihominis RefSeq representative	3133180	GCF_048568405.1	Our current reference
			'accession GCF_019048245.1'  //E. hormaechei RefSeq representative	158836	GCF_019048245.1	Our current reference
		)
		| NCBI_DATASET_DOWNLOAD_GENOME
		
		def db_ch = genomes_ch
			.collect()
			| STEP1_ACCESSIONS_TSV
			| STEP2_MOLECULE_TYPES_TSV

		STEP3_DB_BUILD(
			STEP2_MOLECULE_TYPES_TSV.out,
			NCBI_TAXDUMP_DOWNLOAD.out,
			file("${moduleDir}/assets/db_build.R")
		)

	emit:
		dir   = STEP3_DB_BUILD.out.dir
		fasta = STEP3_DB_BUILD.out.fasta
		tsv   = STEP3_DB_BUILD.out.tsv
}

