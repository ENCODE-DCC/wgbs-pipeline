# Contributing

To work on this pipeline, you will want to have the following installed:
* Docker Desktop
* Java > 1.8
* Python 3.7.4 or higher, tests might fail with lower versions.
  * Install the following:

    ```bash
    $ pip install caper tox
    ```

## Developer Guidelines

Before making a PR, make sure that:
1. `tox` runs with no errors. Invoke it with no arguments like this to perform all checks:
```bash
$ tox
```

2. Make sure that you commit and push any formatting changes that may have occurred above.

## Useful Tips

* To build Docker images manually, run `docker build . -f Dockerfile -t MY_REPO:MY_TAG`
* You can run individual WDL tests manually with `tests/wdl/test.sh tests/wdl/test_task/test_TASK_NAME.wdl tests/wdl/test_task/test_TASK_NAME_input.json. This will invoke `caper` for you. Note you will need to first update the input JSON with `docker` input.
