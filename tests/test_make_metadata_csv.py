from contextlib import contextmanager

import pytest

from wgbs_pipeline.make_metadata_csv import get_parser, main


@contextmanager
def does_not_raise():
    yield


@pytest.mark.parametrize(
    "args,condition,num_files",
    [
        (['-n', 'my_name', '--files', 'foo.fastq.gz'], does_not_raise(), 1),
        (['-n', 'my_name', '--files', 'foo.fastq.gz', "bar.fastq.gz"], does_not_raise(), 2),
        (['--files', 'foo.fastq.gz'], pytest.raises(), 1),
        (['-n', 'my_name'], pytest.raises(), 0),
    ],
)
def test_parser(args, condition, num_files):
    parser = get_parser()
    with condition:
        parsed_args = parser.parse_args(args)
    assert len(parsed_args.files) == num_files
