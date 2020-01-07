import builtins
from contextlib import contextmanager
from typing import List

import attr
import pytest

from wgbs_pipeline.make_metadata_csv import get_parser, main, process


@contextmanager
def does_not_raise():
    yield


@attr.s(auto_attribs=True)
class StubArgs:
    sample_name: str
    files: List[str]


@pytest.mark.parametrize(
    "args,condition,num_files",
    [
        (["-n", "my_name", "--files", "foo.fastq.gz"], does_not_raise(), 1),
        (
            ["-n", "my_name", "--files", "foo.fastq.gz", "bar.fastq.gz"],
            does_not_raise(),
            2,
        ),
        (["--files", "foo.fastq.gz"], pytest.raises(SystemExit), 1),
        (["-n", "my_name"], pytest.raises(SystemExit), 0),
    ],
)
def test_parser(args, condition, num_files):
    parser = get_parser()
    with condition:
        parsed_args = parser.parse_args(args)
        assert len(parsed_args.files) == num_files


@pytest.mark.parametrize(
    "args,condition,expected_row_len",
    [
        (StubArgs("foo", ["file1"]), does_not_raise(), 4),
        (StubArgs("bar", ["file1", "file2"]), does_not_raise(), 5),
        (StubArgs("baz", ["file1", "file2", "file3"]), pytest.raises(ValueError), 0),
    ],
)
def test_process(args, condition, expected_row_len):
    with condition:
        output = process(args)
        assert all(len(row) == expected_row_len for row in output)


def test_main(mocker):
    """
    The assert looks a little wonky here. The first index extracts the args of the call,
    and the second extracts the first positional arg.
    """
    mocker.patch("builtins.open", mocker.mock_open())
    testargs = ["prog", "-n", "foo", "--files", "a.fastq.gz"]
    mocker.patch("sys.argv", testargs)
    main()
    assert builtins.open.call_args[0][0] == "foo_metadata.csv"
