# Contributing

To work on this pipeline, you will want to have the following installed:
* Docker Desktop
* Java > 1.8
* Python 3.7
  * Install the following:

    ```bash
    $ pip install caper tox
    ```
* Rust 1.40.0
  * The easiest way to get up and running is with [rustup](https://rustup.rs/)
  * Once `rustup` is installed, install `rustfmt` and `clippy`:

    ```bash
    $ rustup component add clippy
    $ rustup component add rustfmt
    ```
* R 3.6.2, available from https://cran.r-project.org/
  * Once R is installed, run it and install dev dependencies:

    ```R
    $ install.packages(c("docopt", "styler", "git2r", "lintr")))
    ```

## Developer Guidelines

Before making a PR, make sure that:
1. `tox` runs with no errors. Invoke it with no arguments like this to perform all checks:
```bash
$ tox
```

2. Make sure that you commit and push any formatting changes that may have occurred above.
3. Run `cargo test` and make sure the Rust tests pass.

## Useful Tips

* Running `cargo build` will build all of the Rust binaries if you need to quickly test them. Use `cargo build --release` if you need them to be fast.
* To build Docker images manually, run `docker build . -f Dockerfile -t MY_REPO:MY_TAG`
* You can run individual WDL tests manually with `tests/wdl/test.sh tests/wdl/test_task/test_TASK_NAME.wdl tests/wdl/test_task/test_TASK_NAME_input.json DOCKER_IMAGE`. This will invoke `caper` for you.
