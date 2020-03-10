FROM ubuntu:16.04@sha256:181800dada370557133a502977d0e3f7abda0c25b9bbb035f199f5eb6082a114 as main

LABEL maintainer "Paul Sud"
LABEL maintainer.email "encode-help@lists.stanford.edu"

# Need to add apt key to read to install R 3.6.2: https://stackoverflow.com/a/56378217
# hash -r gets python3 symlink to work
RUN apt-get update && \
    apt-get install -y software-properties-common apt-transport-https && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/" && \
    add-apt-repository -y ppa:deadsnakes/ppa && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    apt-get update && \
    apt-get install -y \
    autoconf \
    build-essential \
    git \
    gsl-bin \
    libbz2-dev \
    libcurl4-openssl-dev \
    libgsl0-dev \
    liblzma-dev \
    libmysqlclient-dev \
    libncurses5-dev \
    libpng-dev \
    libssl-dev \
    libxml2-dev \
    pigz \
    python3.7 \
    r-base=3.6.2-1xenial \
    uuid-dev \
    wget \
    zlib1g-dev \
    && ln -s /usr/bin/python3.7 /usr/local/bin/python3 && hash -r \
    && wget https://bootstrap.pypa.io/get-pip.py && python3 get-pip.py && rm get-pip.py \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir beautifulsoup4 matplotlib==3.0.2 multiprocess

# Install bsmooth.R dependencies
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'https://cloud.r-project.org'; options(repos = r);" > ~/.Rprofile && \
    Rscript -e 'install.packages("BiocManager")' && \
    Rscript -e 'BiocManager::install("bsseq")' && \
    Rscript -e 'install.packages("optparse")'

RUN mkdir /software
WORKDIR /software
ENV PATH="/software:${PATH}"

# Install gemBS
RUN git clone --depth 10 --recursive https://github.com/heathsc/gemBS.git && \
    cd gemBS && git checkout 0724ced33f122376d22c6ab9ef6d5aad548f5595 && \
    python3 setup.py install --user
ENV PATH="/root/.local/bin:${PATH}"

# Instal UCSC v377 bedToBigBed util
RUN git clone https://github.com/ENCODE-DCC/kentUtils_bin_v377 && \
    rm $(find kentUtils_bin_v377/bin/ -type f -not -path '*bedToBigBed')
ENV PATH="${PATH}:/software/kentUtils_bin_v377/bin"
ENV LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/software/kentUtils_bin_v377/lib"

# Compile Rust binaries as a separate stage so we don't bloat pipeline image with Rust
# build toolchain. rustc default target in the image is x86_64-unknown-linux-gnu, which
# is what we want
FROM rust:1.40.0-slim-stretch@sha256:5b30cf90cc7cd948207850d470e56c6bbc3014445de97716af479dd6fc69a763 as builder
WORKDIR /wgbs-pipeline
COPY . .
RUN cargo build --release

FROM main

# Add source files (Python, R)
COPY wgbs_pipeline/*.* wgbs_pipeline/

# Add conf files
COPY conf/* conf/

# Add compiled Rust binaries from other stage of build
COPY --from=builder /wgbs-pipeline/target/release/gembs-to-bismark-bed-converter wgbs_pipeline/
COPY --from=builder /wgbs-pipeline/target/release/bismark-bsmooth-to-encode-bed-converter wgbs_pipeline/

# Add to path
ENV PATH="/software/wgbs_pipeline:${PATH}"
