import builtins
from contextlib import contextmanager
from typing import List

import attr
import pytest

from wgbs_pipeline.make_metadata_csv import get_parser, main, process, read_file_tsv


@contextmanager
def does_not_raise():
    yield


@attr.s(auto_attribs=True)
class StubArgs:
    sample_names: List[str]
    files: str
    barcode_prefix: str = "sample_"


@pytest.mark.parametrize(
    "args,condition,num_names",
    [
        (["-n", "my_name", "--files", "foo.fastq.gz"], does_not_raise(), 1),
        (["-n", "my_name", "my_name2", "--files", "foo.fastq.gz"], does_not_raise(), 2),
        (["--files", "foo.fastq.gz"], pytest.raises(SystemExit), 1),
        (["-n", "my_name"], pytest.raises(SystemExit), 0),
    ],
)
def test_parser(args, condition, num_names):
    parser = get_parser()
    with condition:
        parsed_args = parser.parse_args(args)
        assert len(parsed_args.sample_names) == num_names


@pytest.mark.parametrize(
    "args,parsed_tsv,condition,expected_row_len",
    [
        (StubArgs(["foo"], "path.tsv"), [["f1.fastq.gz"]], does_not_raise(), 4),
        (
            StubArgs(["bar"], "path.tsv"),
            [["f1.fastq.gz", "f2.fastq.gz"]],
            does_not_raise(),
            5,
        ),
        (
            StubArgs(["baz", "qux"], "path.tsv"),
            [["f1.fastq.gz"]],
            pytest.raises(ValueError),
            0,
        ),
    ],
)
def test_process(mocker, args, parsed_tsv, condition, expected_row_len):
    mocker.patch(
        "wgbs_pipeline.make_metadata_csv.read_file_tsv", return_value=parsed_tsv
    )
    with condition:
        output = process(args)
        assert all(len(row) == expected_row_len for row in output)


@pytest.mark.parametrize(
    "tsv_contents,condition,expected",
    [
        (
            "/a/f1.fastq.gz\t/b/f2.fastq.gz\n",
            does_not_raise(),
            [["f1.fastq.gz", "f2.fastq.gz"]],
        ),
        ("f1.fastq.gz\tf2.fastq.gz\tf3.fastq.gz\n", pytest.raises(ValueError), [[]]),
    ],
)
def test_read_file_tsv(mocker, tsv_contents, condition, expected):
    """
    mock_open does not support iteration (fixed in Python 3.8 but not backported), so we
    need to patch it. See https://bugs.python.org/issue21258 . We can resolve this by
    supplying a __iter__ method to the mock: https://stackoverflow.com/a/24779923
    """
    mocker.patch("builtins.open", mocker.mock_open(read_data=tsv_contents))
    mocker.patch(
        "builtins.open.return_value.__iter__", lambda self: iter(self.readline, "")
    )
    with condition:
        result = read_file_tsv("foo")
        assert result == expected


def test_main(mocker):
    """
    The assert looks a little wonky here. The first index extracts the args of the call,
    and the second extracts the first positional arg. We use the same trick of patching
    open()'s __iter__ so that the csv library can loop over it.
    """
    mocker.patch("builtins.open", mocker.mock_open(read_data="f1.fastq.gz"))
    mocker.patch(
        "builtins.open.return_value.__iter__", lambda self: iter(self.readline, "")
    )
    testargs = ["prog", "-n", "foo", "--files", "path.tsv"]
    mocker.patch("sys.argv", testargs)
    main()
    assert builtins.open.call_args[0][0] == "metadata.csv"
