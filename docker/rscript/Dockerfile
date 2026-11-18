FROM docker.io/bioconductor/bioconductor_docker:RELEASE_3_23-R-4.6.0

RUN apt-get update && apt-get install -y \
      jq \
    && rm -rf /var/cache/apt/* /var/lib/apt/lists/*;

RUN install2.r --error --skipinstalled --ncpus -1 \
        optparse tidyverse \
        tidygraph ggraph openxlsx patchwork \
        rmarkdown gt flexdashboard plotly reactable visNetwork \
        jqr \
    && rm -rf /tmp/downloaded_packages
RUN Rscript -e 'BiocManager::install(c("Biostrings","GenomicRanges","BiocParallel","GenomicAlignments","rtracklayer"))'

VOLUME /cwd
WORKDIR /cwd
ENTRYPOINT []
