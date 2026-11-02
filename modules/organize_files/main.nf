nextflow.enable.types = true

process ORGANIZE_FILES {
    memory 2.GB
    cpus 1
    time 1.h
    executor 'local'

    input:
    tuple(meta: Object, file_pairs: List) // list of [src, target]; src is a Path (symlinked) or a value (written as text)

    output:
    tuple(meta, file('output'))

    script:
    def cmd = file_pairs
        .findAll({ p -> p[0] })
        .collect({ p ->
            def src = p[0]
            def dest_path = "output/${p[1]}"
            src instanceof Path
                ? "mkdir -p \$(dirname '${dest_path}') && ln -s '${src}' '${dest_path}'"
                : "mkdir -p \$(dirname '${dest_path}') && echo '${src}' > '${dest_path}'"
        })
    """
    trap '' PIPE
    mkdir -p output
    ${cmd.join('\n')}
    """
}
