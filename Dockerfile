FROM ubuntu:16.04@sha256:181800dada370557133a502977d0e3f7abda0c25b9bbb035f199f5eb6082a114

LABEL maintainer "Paul Sud"
LABEL maintainer.email "encode-help@lists.stanford.edu"

# hash -r gets python3 symlink to work
RUN apt-get update && \
    apt-get install -y software-properties-common apt-transport-https && \
    add-apt-repository -y ppa:deadsnakes/ppa && \
    apt-get update && \
    apt-get install -y \
    autoconf \
    build-essential \
    git \
    gsl-bin \
    jq \
    lbzip2 \
    libbz2-dev \
    libcurl4-openssl-dev \
    libgsl0-dev \
    liblzma-dev \
    libncurses5-dev \
    libssl-dev \
    libxml2-dev \
    python3.7 \
    wget \
    zlib1g-dev \
    && ln -s /usr/bin/python3.7 /usr/local/bin/python3 && hash -r \
    && wget https://bootstrap.pypa.io/get-pip.py && python3 get-pip.py && rm get-pip.py \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir /software
WORKDIR /software
ENV PATH="/software:${PATH}"

# Install gemBS
RUN git clone --recursive https://github.com/heathsc/gemBS.git && \
    cd gemBS && \
    git checkout 3bed486a707671e452642b4e357d9316e36db22b && \
    python3 setup.py install --user && \
    cd .. && \
    rm -rf gemBS /usr/local/build
ENV PATH="/root/.local/bin:${PATH}"

# Install UCSC v377 bedToBigBed util
RUN git clone https://github.com/ENCODE-DCC/kentUtils_bin_v377 && \
    rm $(find kentUtils_bin_v377/bin/ -type f -not -path '*bedToBigBed' -not -path '*bedGraphToBigWig')
ENV PATH="${PATH}:/software/kentUtils_bin_v377/bin"
ENV LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/software/kentUtils_bin_v377/lib"

# Install xsv for selecting out bed columns for generating coverage bigiwg
RUN wget https://github.com/BurntSushi/xsv/releases/download/0.13.0/xsv-0.13.0-x86_64-unknown-linux-musl.tar.gz && \
    tar xf xsv-0.13.0-x86_64-unknown-linux-musl.tar.gz && \
    mv xsv /usr/local/bin && \
    rm xsv-0.13.0-x86_64-unknown-linux-musl.tar.gz

COPY requirements*.txt ./
RUN pip3 install --no-cache-dir -r requirements.txt -r requirements-gembs.txt && \
    rm requirements*.txt

# Add source files (Python, R)
COPY wgbs_pipeline/*.* wgbs_pipeline/

# Add conf files
COPY conf/* conf/

# Add to path
ENV PATH="/software/wgbs_pipeline:${PATH}"

# pytest-workflow tests don't run as root, need to chmod dependencies
RUN chmod -R a+rwx /software /root
