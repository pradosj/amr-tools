FROM docker.io/python:3.10-slim-bullseye

RUN apt-get update -qq && apt-get install -y -qq \
    git \
    build-essential \
    libz-dev \
    ncbi-blast+ \
    curl \
    && rm -rf /var/cache/apt/* /var/lib/apt/lists/*

# KMA
RUN git -C /usr/src/ clone --depth 1 https://bitbucket.org/genomicepidemiology/kma.git \
    && cd /usr/src/kma \
    && make \
    && cp /usr/src/kma/kma /usr/src/kma/kma_index /usr/bin

# plasmidfinder
RUN mkdir -p /db \
    && git -C /db clone --depth 1 https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git \
    && cd /db/plasmidfinder_db && python3 INSTALL.py
ENV PLASMID_DB=/db/plasmidfinder_db
RUN git -C /usr/src/ clone --depth 1 https://bitbucket.org/genomicepidemiology/plasmidfinder/ \
    && cp /usr/src/plasmidfinder/plasmidfinder.py /usr/bin


# mlst
RUN mkdir -p /db \
    && git -C /db clone --depth 1 https://bitbucket.org/genomicepidemiology/mlst_db.git \
    && cd /db/mlst_db && python3 INSTALL.py
RUN git -C /usr/src/ clone --depth 1 https://bitbucket.org/genomicepidemiology/mlst/ \
    && cp /usr/src/mlst/mlst.py /usr/bin \
    && pip install cgecore==1.5.5 tabulate==0.9.0 biopython


# resfinder
RUN mkdir -p /db \
    && git -C /db clone --depth 1 https://git@bitbucket.org/genomicepidemiology/resfinder_db.git \
    && git -C /db clone --depth 1 https://git@bitbucket.org/genomicepidemiology/pointfinder_db.git \
    && git -C /db clone --depth 1 https://git@bitbucket.org/genomicepidemiology/disinfinder_db.git \
    && cd /db/resfinder_db && python3 INSTALL.py \
    && cd /db/pointfinder_db && python3 INSTALL.py \
    && cd /db/disinfinder_db && python3 INSTALL.py
ENV CGE_RESFINDER_RESGENE_DB=/db/resfinder_db/
ENV CGE_RESFINDER_RESPOINT_DB=/db/pointfinder_db/
ENV CGE_DISINFINDER_DB=/db/disinfinder_db/
RUN curl -sSL https://pdm-project.org/install-pdm.py | python3 - --path /usr/local \
    && git -C /usr/src/ clone --depth 1 https://bitbucket.org/genomicepidemiology/resfinder.git \
    && cd /usr/src/resfinder && pdm build && pip3 install dist/*.whl



