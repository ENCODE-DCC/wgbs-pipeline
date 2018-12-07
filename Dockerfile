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
	jq

RUN pip3 install matplotlib multiprocess

# Make directory for all  softwares
RUN mkdir /software
WORKDIR /software
ENV PATH="/software:${PATH}"


# Install gemBS
RUN git clone --recursive https://github.com/heathsc/gemBS.git
RUN cd gemBS && python3 setup.py install --user
ENV PATH="/root/.local/bin:${PATH}"

# Set up Python helpers for WDL scripts
RUN mkdir -p helpers
COPY helpers helpers
RUN chmod +x -R helpers/ 
ENV PATH="/software/helpers:${PATH}"

ENTRYPOINT ["/bin/bash","-c"]