

docker build --platform linux/amd64 -t registry.gitlab.unige.ch/amr-genomics/rscript:v3 ./
docker push registry.gitlab.unige.ch/amr-genomics/rscript:v3

# Tests
docker run --rm -it registry.gitlab.unige.ch/amr-genomics/rscript:v3 bash




