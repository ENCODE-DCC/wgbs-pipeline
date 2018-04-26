FROM ubuntu:16.04

LABEL maintainer="Ulugbek Baymuradov"

RUN apt-get update && apt-get install -y \
	apt-utils \
	libc-dev \
	build-essential \
	software-properties-common \
	libncurses5-dev \
    libncursesw5-dev \
    libfreetype6-dev \
    libgsl0-dev \
    r-base \
    pigz \
	python-pip \
	python3-dev \
    python3-pip \
	git \
	wget \
	unzip \
	p7zip-full \
	pkg-config

# Make directory for all softwares
RUN mkdir /software
WORKDIR /software
ENV PATH="/software:${PATH}"

# Install OpenBLAS 0.2.19 (NEHALEM architecture, single-threaded)
# to have the same result as of Conda's non-MKL openblas 0.2.19
# setting TARGET=HASWELL will give you different result.
# so let's stick to TARGET=NEHALEM
# RUN git clone --branch v0.2.19 https://github.com/xianyi/OpenBLAS
# RUN cd OpenBLAS && make FC=gfortran TARGET=NEHALEM USE_THREAD=0 && make PREFIX=/opt/openblas install
# ENV LD_LIBRARY_PATH="/opt/openblas/lib:${LD_LIBRARY_PATH}"


# Install samtools dependencies
RUN wget http://zlib.net/zlib-1.2.11.tar.gz && tar -xvf zlib-1.2.11.tar.gz
RUN cd zlib-1.2.11 && ./configure && make && make install
RUN wget http://bzip.org/1.0.6/bzip2-1.0.6.tar.gz && tar -xvf bzip2-1.0.6.tar.gz
RUN cd bzip2-1.0.6 && make && make install
RUN wget https://tukaani.org/xz/xz-5.2.3.tar.gz && tar -xvf xz-5.2.3.tar.gz
RUN cd xz-5.2.3 && ./configure && make && make install

# Upgrade pip
RUN pip install --upgrade pip

# Install setuptools
RUN pip install -U setuptools


# Install basic python2 packages
RUN pip install common python-dateutil cython
RUN pip3 install common python-dateutil cython

# Install numpy 1.11.3 (python2/3)
# RUN git clone --branch v1.11.3 https://github.com/numpy/numpy
# COPY /conf/site.cfg numpy/
# RUN cd numpy && python setup.py install

# Install numpy via pip
RUN pip install numpy

# Install samtools 1.3
RUN git clone --branch 1.3 https://github.com/samtools/samtools.git
RUN git clone --branch 1.3 git://github.com/samtools/htslib.git
RUN cd samtools && make && make install

# Install bcftools 1.3.1
RUN git clone --branch 1.3.1 git://github.com/samtools/bcftools.git
RUN cd bcftools && make

# Install matplotlib
RUN pip install --upgrade matplotlib[mplot3d]

# Install Texlive, this takes long time and adds 10+ GB to image size
# We are going to comment it out in development

# RUN wget http://mirror.hmc.edu/ctan/systems/texlive/Images/texlive.iso
# RUN 7z x -oextracted_texlive texlive.iso
# RUN cd extracted_texlive && chmod +x install-tl && echo 'I' | ./install-tl
# RUN PATH="/usr/local/texlive/2017/bin/x86_64-linux:${PATH}"

# Get wigToBigWig
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
RUN chmod +x wigToBigWig

# Install gemBS
RUN git clone --recursive https://github.com/heathsc/gemBS.git
RUN sed -i 's+^GSL_LIB =.*+GSL_LIB = -L/usr/lib/x86_64-linux-gnu/+g' /software/gemBS/tools/bs_call/Gsl.mk
RUN sed -i 's+^GSL_INC =.*+GSL_INC = -I/usr/include/gsl/+g' /software/gemBS/tools/bs_call/Gsl.mk
RUN cd gemBS && python setup.py install --user
ENV PATH="/root/.local/bin:${PATH}"

# Set up test data
RUN cd gemBS/test && tar -xvf example.tar.gz 

# Set up Python helpers for WDL scripts
RUN mkdir helpers/
COPY helpers/ helpers/
ENV PATH="/software/helpers:${PATH}"

ENTRYPOINT ["/bin/bash","-c"]