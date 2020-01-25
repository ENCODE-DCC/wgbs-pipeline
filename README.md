# wgbs-pipeline

[![CircleCI](https://circleci.com/gh/ENCODE-DCC/wgbs-pipeline.svg?style=svg)](https://circleci.com/gh/ENCODE-DCC/wgbs-pipeline)
[![Docker Repository on Quay](https://quay.io/repository/encode-dcc/wgbs-pipeline/status "Docker Repository on Quay")](https://quay.io/repository/encode-dcc/wgbs-pipeline)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![MIT License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

## Overview

An [ENCODE](https://www.encodeproject.org/) pipeline for processing whole-genome bisulfite sequencing (WGBS) data using [gemBS](https://github.com/heathsc/gemBS) for alignment and methylation extraction and [Bsmooth](https://bioconductor.org/packages/release/bioc/html/bsseq.html) for smoothing of CpG methylation values.

This pipeline is currently in prerelease, and as such may be unstable.

## Installation

1) Git clone this pipeline.
    ```bash
    $ git clone https://github.com/ENCODE-DCC/wgbs-pipeline
    ```

2) Install [Caper](https://github.com/ENCODE-DCC/caper), requires `java` > 1.8 and `python` > 3.4.1 . Caper is a python wrapper for [Cromwell](https://github.com/broadinstitute/cromwell).
    ```bash
    $ pip install caper  # use pip3 if it doesn't work
    ```

3) Follow [Caper's README](https://github.com/ENCODE-DCC/caper) carefully to configure it for your platform (local, cloud, cluster, etc.)
    > **IMPORTANT**: Configure your Caper configuration file `~/.caper/default.conf` correctly for your platform.

## Usage

```bash
$ caper run wgbs-pipeline.wdl -i input.json -m testrun_metadata.json
```

## Description of files

The gemBS index file is a .tar.gz archive containing the following files produced by running gemBS index:
* ${prefix}.gemBS.contig_md5
* ${prefix}.gemBS.ref.gzi
* ${prefix}.BS.gem
* ${prefix}.gemBS.ref
* ${prefix}.BS.info
* ${prefix}.contig.sizes
* ${prefix}.gemBS.ref.fai

## Contributing

We welcome comments, questions, suggestions, bug reports, feature requests, and pull requests (PRs). Please use one of the existing Github issue templates if applicable. When contributing code, please follow the [Developer Guidelines](docs/CONTRIBUTING.md#developer-guidelines).
