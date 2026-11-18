

# Introduction

`cgetools` is s docker container with several tools and databases from CGE installed inside.


# Usage

To use the tool, you need to mount the working directory on `/cwd`. Example:
```bash
docker run --rm -v .:/cwd registry.gitlab.unige.ch/amr-genomics/cgetools:main plasmidfinder.py
docker run --rm -v .:/cwd registry.gitlab.unige.ch/amr-genomics/cgetools:main mlst.py
docker run --rm -v .:/cwd registry.gitlab.unige.ch/amr-genomics/cgetools:main python3 -m resfinder
```


# Building

To use the tool, you need to mount the working directory on `/cwd`. Example:
```bash
docker build ./ --platform linux/amd64 -t registry.gitlab.unige.ch/amr-genomics/cgetools:v2.0
docker push registry.gitlab.unige.ch/amr-genomics/cgetools:v2.0
```
