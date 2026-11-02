nextflow.enable.types = true

process RMD_RENDER {
    container "registry.gitlab.unige.ch/amr-genomics/rscript:v2"
    label 'rscript'
    memory 8.GB
    cpus 2
    time 30.m

    input:
    // inputs: files and folders used by the template, staged under their own name
    // render_params: R expressions of the render params, e.g. "indir='folder'"
    record(template: Path, inputs: List<Path>, render_params: String?)

    output:
    file('report.html')

    script:
    """
    #!/usr/bin/env Rscript
    rmarkdown::render(
        knit_root_dir = getwd(),
        '${template}',
        ${render_params ? "params = list(${render_params})," : ''}
        output_dir = getwd(),
        output_file = "report.html"
    )
    """

    stub:
    """
    touch report.html
    """
}
