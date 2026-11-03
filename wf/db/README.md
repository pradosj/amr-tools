# db

`db` a nextflow pipeline to build the databases of the `amr-tools` pipelines in a folder, used as a cache with `--db`:

- `checkm2/`: CheckM2 database (assembly completeness, `ASSEMBLIES_QC`);
- `org_db/`: orgfinder reference genomes, downloaded from NCBI (organism detection, `AMR_ANNOT`);
- `amrfinderplus/`: AMRFinderPlus database (`AMR_ANNOT`).

It also makes a tar archive of the folder (`amr-tools-db.tar.gz`), that the pipelines take with `--db` as a path or an URL.

# Usage

```bash
# build the databases in /data/amr-tools-db (folders and archive)
nextflow run pradosj/amr-tools/wf/db -output-dir /data/amr-tools-db

# use them
nextflow run pradosj/amr-tools --db=/data/amr-tools-db --reads.ont='ont/*.fastq.gz'
nextflow run pradosj/amr-tools --db=/data/amr-tools-db/amr-tools-db.tar.gz --reads.ont='ont/*.fastq.gz'
```

Options:
- `--db`: an existing databases folder (or archive) to complete; only the missing databases are built. With `--db` and `-output-dir` set to the same folder, the folder is completed in place.
- `--archive`: name of the archive (`<archive>.tar.gz`, default `amr-tools-db`); `--archive=''` for no archive. To update a database, delete its folder and run the pipeline again.

# Test the pipeline

```bash
nextflow run pradosj/amr-tools/wf/db/main.nf -profile test -stub
```
