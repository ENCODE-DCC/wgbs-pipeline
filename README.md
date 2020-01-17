# wgbs-pipeline

An [ENCODE](https://www.encodeproject.org/) pipeline for processing whole-genome bisulfite sequencing (WGBS) data using [gemBS](https://github.com/heathsc/gemBS) pipeline for alignment and methylation extraction and [Bsmooth](https://bioconductor.org/packages/release/bioc/html/bsseq.html) for smoothing of CpG methylation values.

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

## Contributing

We welcome comments, questions, suggestions, bug reports, feature requests, and pull requests (PRs). Before making PRs, please follow the [PR Guidelines](https://github.com/ENCODE-DCC/wgbs-pipeline#pr-guidelines).

## PR Guidelines

Before making a PR, make sure that:
1. `tox` passes. To run, ensure you have `python` 3.7 available on your `$PATH`, install `tox` with `pip install tox`, and run `tox` (no arguments) from the repo root.
2. (only if Rust code changes) Lint and format Rust files by running the following commands (you can set these all to run automatically in an editor like VS Code). To run, first install [rustup](https://rustup.rs/), which will install `cargo`, then install `rustfmt` (`rustup component add rustfmt`) and `clippy` (`rustup component add clippy`)
    ```bash
    $ cargo fmt
    $ cargo check
    $ cargo clippy
    ```
3. (only if R code changes) Format R files by loading them in [RStudio free desktop version](https://rstudio.com/products/rstudio/download/#download) then clicking `Code` > `Reformat Code` or pressing `Shift` + `Cmd` + `A` (Mac) and saving.
4. Make sure that you commit and push any formatting changes after any formatters.
