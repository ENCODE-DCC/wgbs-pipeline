FROM ubuntu:16.04

LABEL maintainer="Ulugbek Baymuradov"

RUN apt-get update && apt-get install -y \
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
	libpng-dev \
	uuid-dev \
	libmysqlclient-dev \
	jq

RUN pip3 install matplotlib==3.0.2 multiprocess

# Make directory for all  softwares
RUN mkdir /software
WORKDIR /software
ENV PATH="/software:${PATH}"


# Install gemBS
# We do this dummy add step to invalidate the cache if the master HEAD ref changes
# Otherwise the we won't pull the latest code changes
# See https://stackoverflow.com/a/39278224
ADD https://api.github.com/repos/heathsc/gemBS/git/refs/heads/master versions.json
RUN git clone --recursive https://github.com/heathsc/gemBS.git
RUN cd gemBS && python3 setup.py install --user
ENV PATH="/root/.local/bin:${PATH}"

COPY wgbs_pipeline wgbs_pipeline
ENV PATH="/software/wgbs_pipeline:${PATH}"

ENTRYPOINT ["/bin/bash","-c"]
