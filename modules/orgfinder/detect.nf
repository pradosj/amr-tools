nextflow.enable.types = true

process ORGFINDER_DETECT {
    container 'registry.gitlab.unige.ch/amr-genomics/orgfinder:v0.2'
    memory 20.GB
    cpus 8
    time 15.m

    input:
    record(meta: Map, fasta: Path, args: String?)

    output:
    record(
        meta: meta,
        dir: file('orgfinder'),
        org_name: env('ORG_NAME'),
        org_acc: env('ORG_ACC'),
        ani: env('ORG_ANI')
    )

    script:
    """
    mkdir ./orgfinder && cp /app/db/tax.tsv /app/db/tax.rds /app/db/db.tsv ./orgfinder
    orgfinder ${args ?: ''} --thread=${task.cpus} --out 'orgfinder/ani.tsv' '${fasta}'
    Rscript -e 'readr::read_tsv("orgfinder/ani.tsv",show_col_types=FALSE) |>
        dplyr::slice_max(ANI,n=1,with_ties = FALSE) |>
        dplyr::mutate(stringr::str_glue("ORG_ACC=\\"{assembly_acc}\\"\\nORG_NAME=\\"{org_name}\\"\\nORG_ANI=\\"{ANI}\\"")) |>
        dplyr::pull() |>
        writeLines("org.env")
    '
    source "org.env"
    """

    stub:
    """
    mkdir -p orgfinder
    touch orgfinder/ani.tsv orgfinder/tax.tsv orgfinder/db.tsv
    export ORG_NAME="xxx"
    export ORG_ACC="xxx"
    export ORG_ANI="90"
    """
}
