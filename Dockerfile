FROM ubuntu:16.04@sha256:181800dada370557133a502977d0e3f7abda0c25b9bbb035f199f5eb6082a114 as main

LABEL maintainer "Paul Sud"
LABEL maintainer.email "encode-help@lists.stanford.edu"

# Need to add apt key to read to install R 3.6.2: https://stackoverflow.com/a/56378217
RUN apt-get update && \
    apt-get install -y software-properties-common apt-transport-https && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/" && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    apt-get update && \
    apt-get install -y \
    python3 \
    build-essential \
    git \
    python3-pip \
    wget \
    pigz \
    zlib1g-dev \
    libbz2-dev \
    gsl-bin \
    libgsl0-dev \
    libncurses5-dev \
    liblzma-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libpng-dev \
    uuid-dev \
    libmysqlclient-dev \
    r-base=3.6.2-1xenial \
    jq \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir matplotlib==3.0.2 multiprocess

# Install bsmooth.R dependencies
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'https://cloud.r-project.org'; options(repos = r);" > ~/.Rprofile && \
    Rscript -e 'install.packages("BiocManager")' && \
    Rscript -e 'BiocManager::install("bsseq")' && \
    Rscript -e 'install.packages("optparse")'

RUN mkdir /software
WORKDIR /software
ENV PATH="/software:${PATH}"

# Compile Rust binaries as a separate stage so we don't bloat pipeline image with Rust
# build toolchain. rustc default target in the image is x86_64-unknown-linux-gnu, which
# is what we want
FROM rust:1.40.0-slim-stretch as builder
WORKDIR /wgbs-pipeline
COPY . .
RUN cargo build --release

FROM main

# Install gemBS
RUN git clone --depth 1 --recursive https://github.com/heathsc/gemBS.git && \
    cd gemBS && git checkout 6d7a8ab25c2c44e6c6cca1485bba5b5fbaafc88f && \
    python3 setup.py install --user
ENV PATH="/root/.local/bin:${PATH}"

# Instal Kent bedToBigBed util
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed && \
    chmod +x bedToBigBed

# Add source files (Python, R)
COPY wgbs_pipeline/*.* wgbs_pipeline/

# Add compiled Rust binaries from other stage of build
COPY --from=builder /wgbs-pipeline/target/release/gembs-to-bismark-bed-converter wgbs_pipeline/
COPY --from=builder /wgbs-pipeline/target/release/bismark-bsmooth-to-encode-bed-converter wgbs_pipeline/

# Add to path
ENV PATH="/software/wgbs_pipeline:${PATH}"
