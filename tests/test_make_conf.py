import builtins
from contextlib import suppress as does_not_raise
from typing import List, Optional, Tuple

import attr
import pytest

from wgbs_pipeline.make_conf import get_parser, main, make_conf, write_conf


@attr.s(auto_attribs=True)
class StubArgs:
    reference: str
    extra_reference: str
    underconversion_sequence: Optional[str] = None
    include_file: Optional[str] = None
    outfile: str = "gembs.conf"
    num_threads: int = 8
    num_jobs: int = 3


@pytest.mark.parametrize(
    "args,condition",
    [
        (
            [
                "-t",
                "3",
                "-j",
                "8",
                "-r",
                "ref.fa.gz",
                "-e",
                "extra_ref.fa.gz",
                "-u",
                "chrL",
                "-i",
                "include.conf",
                "-o",
                "my.conf",
            ],
            does_not_raise(),
        ),
        (["-r", "ref.fa.gz", "-e", "extra_ref.fa.gz"], does_not_raise()),
        (["-r", "ref.fa.gz"], pytest.raises(SystemExit)),
        (["-e", "extra.fa.gz"], pytest.raises(SystemExit)),
    ],
)
def test_parser(args: List[str], condition):
    parser = get_parser()
    with condition:
        parser.parse_args(args)


@pytest.mark.parametrize(
    "args,expected_num_rows,assertions",
    [
        (
            StubArgs("/foo/reference.fa.gz", "/bar/exref.fa.gz"),
            11,
            [
                (0, "reference = reference/reference.fa.gz"),
                (1, "extra_references = reference/exref.fa.gz"),
                (4, "sequence_dir = ${base}/fastq/@SAMPLE"),
                (9, "threads = 8"),
                (10, "jobs = 3"),
            ],
        ),
        (
            StubArgs("reference.fa.gz", "exref.fa.gz", underconversion_sequence="chrL"),
            13,
            [(11, "[mapping]"), (12, "underconversion_sequence = chrL")],
        ),
        (
            StubArgs(
                "reference.fa.gz", "exref.fa.gz", include_file="my_dir/include.conf"
            ),
            13,
            [(11, "[mapping]"), (12, "include my_dir/include.conf")],
        ),
        (
            StubArgs(
                "/foo/reference.fa.gz", "/bar/exref.fa.gz", num_threads=1, num_jobs=1
            ),
            11,
            [(9, "threads = 1"), (10, "jobs = 1")],
        ),
        (
            StubArgs(
                "reference.fa.gz",
                "exref.fa.gz",
                underconversion_sequence="chrL",
                include_file="include.conf",
            ),
            14,
            [
                (11, "[mapping]"),
                (12, "underconversion_sequence = chrL"),
                (13, "include include.conf"),
            ],
        ),
    ],
)
def test_make_conf(args, expected_num_rows: int, assertions: List[Tuple[int, str]]):
    output = make_conf(args)
    assert len(output) == expected_num_rows
    for row_index, expected_value in assertions:
        assert output[row_index] == expected_value


def test_write_conf(mocker):
    """
    The first two calls to mock_open are the arguments to the open() (filename and mode)
    and the context manager __enter__(). The third call then contains the line of data
    written. The positional arguments to this call are in the second entry of the tuple
    (first is the name of the callee). The last index selects the first (and only)
    positional arg.
    """
    mocker.patch("builtins.open", mocker.mock_open())
    write_conf("my_outfile.conf", ["foo = bar"])
    assert builtins.open.mock_calls[2][1][0] == "foo = bar\n"


def test_main(mocker):
    """
    The assert looks a little wonky here. The first index extracts the args of the call,
    and the second extracts the first positional arg.
    """
    mocker.patch("builtins.open", mocker.mock_open())
    testargs = ["prog", "-r", "ref.fa.gz", "-e", "exref.fa.gz", "-o", "conf.conf"]
    mocker.patch("sys.argv", testargs)
    main()
    assert builtins.open.call_args[0][0] == "conf.conf"
