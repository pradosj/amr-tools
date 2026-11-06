nextflow.enable.types = true

// Render a quarto template; `files` and `extra` are staged next to the template, `args` holds quarto options (e.g. -P name:value)
process QUARTO_RENDER {
    container "registry.gitlab.unige.ch/amr-genomics/rscript:v2"
    label 'rscript'
    memory 8.GB
    cpus 2
    time 30.m

    input:
    record(meta: Record, template: Path, extra: Path?, files: List<Path>, args: String?)

    stage:
    stageAs(template, 'report_template.qmd')

    output:
    record(meta: meta, html: file('report_template.html'))

    script:
    """
    HOME=/tmp quarto render 'report_template.qmd' ${args ?: ''}
    """

    stub:
    """
    touch report_template.html
    """
}
