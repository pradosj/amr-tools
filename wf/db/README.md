# db

`db` a nextflow pipeline to build the databases of the `amr-tools` pipelines in their databases folder (`results/db` by default):

- `checkm2/`: CheckM2 database (assembly completeness, `ASSEMBLIES_QC`);
- `org_db/`: orgfinder reference genomes, downloaded from NCBI (organism detection, `AMR_ANNOT`);
- `amrfinderplus/`: AMRFinderPlus database (`AMR_ANNOT`).

It also makes a tar archive of the folder (`amr-tools-db.tar.gz`), that the pipelines take with `--db.basedir` as a path or an URL.

# Usage

```bash
# build the databases in results/db (and results/amr-tools-db.tar.gz), the default databases folder of the pipelines
nextflow run pradosj/amr-tools/wf/db

# build them elsewhere, and use them
nextflow run pradosj/amr-tools/wf/db -output-dir /data/amr-tools
nextflow run pradosj/amr-tools --db.basedir=/data/amr-tools/db --reads.ont='ont/*.fastq.gz'
nextflow run pradosj/amr-tools --db.basedir=/data/amr-tools/amr-tools-db.tar.gz --reads.ont='ont/*.fastq.gz'
```

Options:
- `--db.basedir` (default `results/db`), `--db.checkm2`, `--db.org_db`, `--db.amrfinderplus`: the databases folder and the paths of the databases in it, as in the other pipelines. Only the databases missing from the folder are built, and they are published to `<output dir>/db/` under their path: with the default `--db.basedir` and `-output-dir`, the folder is completed in place.
- `--archive`: name of the archive (`<archive>.tar.gz`, default `amr-tools-db`); `--archive=''` for no archive. To update a database, delete its folder and run the pipeline again.

# Test the pipeline

```bash
nextflow run pradosj/amr-tools/wf/db/main.nf -profile test -stub
```
