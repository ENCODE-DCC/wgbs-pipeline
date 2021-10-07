import builtins
from contextlib import suppress as does_not_raise
from typing import List

import attr
import pytest

from wgbs_pipeline.make_metadata_csv import get_parser, main, process, read_file_json


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
    "args,parsed_json,condition,expected",
    [
        (
            StubArgs(["foo"], "path.json"),
            [[["f1.fastq.gz"]]],
            does_not_raise(),
            [
                ["Barcode", "Name", "Dataset", "File"],
                ["sample_foo", "foo", "0", "f1.fastq.gz"],
            ],
        ),
        (
            StubArgs(["bar"], "path.json"),
            [[["f1.fastq.gz", "f2.fastq.gz"]]],
            does_not_raise(),
            [
                ["Barcode", "Name", "Dataset", "File1", "File2"],
                ["sample_bar", "bar", "0", "f1.fastq.gz", "f2.fastq.gz"],
            ],
        ),
        (
            StubArgs(["baz"], "path.json"),
            [[["f1.fastq.gz", "f2.fastq.gz"], ["f3.fastq.gz", "f4.fastq.gz"]]],
            does_not_raise(),
            [
                ["Barcode", "Name", "Dataset", "File1", "File2"],
                ["sample_baz", "baz", "0", "f1.fastq.gz", "f2.fastq.gz"],
                ["sample_baz", "baz", "1", "f3.fastq.gz", "f4.fastq.gz"],
            ],
        ),
        (
            StubArgs(["baz", "qux"], "path.json"),
            [[["f1.fastq.gz", "f2.fastq.gz"]], [["f3.fastq.gz", "f4.fastq.gz"]]],
            does_not_raise(),
            [
                ["Barcode", "Name", "Dataset", "File1", "File2"],
                ["sample_baz", "baz", "0", "f1.fastq.gz", "f2.fastq.gz"],
                ["sample_qux", "qux", "0", "f3.fastq.gz", "f4.fastq.gz"],
            ],
        ),
        (
            StubArgs(["baz", "qux"], "path.json"),
            [[["f1.fastq.gz"]]],
            pytest.raises(ValueError),
            [[]],
        ),
    ],
)
def test_process(mocker, args, parsed_json, condition, expected):
    mocker.patch(
        "wgbs_pipeline.make_metadata_csv.read_file_json", return_value=parsed_json
    )
    with condition:
        output = process(args)
        assert output == expected


@pytest.mark.parametrize(
    "json_contents,condition,expected",
    [
        (
            '[[["/a/f1.fastq.gz","/b/f2.fastq.gz"]]]',
            does_not_raise(),
            [[["f1.fastq.gz", "f2.fastq.gz"]]],
        ),
        (
            '[[["/a/f1.fastq.gz","/b/f2.fastq.gz"],["/a/f3.fastq.gz","/b/f4.fastq.gz"]]]',
            does_not_raise(),
            [[["f1.fastq.gz", "f2.fastq.gz"], ["f3.fastq.gz", "f4.fastq.gz"]]],
        ),
        (
            '[[["/a/f1.fastq.gz", "/b/f2.fastq.gz"]],[["/a/f3.fastq.gz", "/b/f4.fastq.gz"]]]',
            does_not_raise(),
            [[["f1.fastq.gz", "f2.fastq.gz"]], [["f3.fastq.gz", "f4.fastq.gz"]]],
        ),
        (
            '[[["f1.fastq.gz", "f2.fastq.gz", "f3.fastq.gz"]]]',
            pytest.raises(ValueError),
            [[[]]],
        ),
        (
            '[[["f1.fastq.gz"], ["f2.fastq.gz", "f3.fastq.gz"]]]',
            pytest.raises(ValueError),
            [[[]]],
        ),
    ],
)
def test_read_file_json(mocker, json_contents, condition, expected):
    mocker.patch("builtins.open", mocker.mock_open(read_data=json_contents))
    with condition:
        result = read_file_json("foo")
        assert result == expected


def test_main(mocker):
    """
    The assert looks a little wonky here. The first index extracts the args of the call,
    and the second extracts the first positional arg.
    """
    mocker.patch("builtins.open", mocker.mock_open(read_data=b'["f1.fastq.gz"]'))
    testargs = ["prog", "-n", "foo", "--files", "path.tsv"]
    mocker.patch("sys.argv", testargs)
    main()
    assert builtins.open.call_args[0][0] == "metadata.csv"
