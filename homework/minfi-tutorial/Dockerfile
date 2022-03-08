FROM rocker/rstudio:4.1.2
ARG DEBIAN_FRONTEND="noninteractive"
RUN install2.r BiocManager
RUN apt-get update \
    && apt-get install --yes \
        libbz2-dev \
        liblzma-dev \
        libxml2 \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*
RUN R -e 'BiocManager::install(c("minfi", "minfiData", "sva"))' \
    && rm -rf /tmp/*
