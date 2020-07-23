import pytest
from qc_utils import QCMetricRecord

from wgbs_pipeline.calculate_average_coverage import (
    calculate_average_coverage,
    calculate_genome_size,
    get_samtools_stats,
    make_qc_record,
)


@pytest.mark.filesystem
def test_get_samtools_stats(mocker):
    mocker.patch("pysam.stats", return_value="SN\tfoo:\t3\n")
    result = get_samtools_stats("bam_path", threads=3)
    assert result == {"foo": 3}


def test_calculate_genome_size(mocker):
    mocker.patch("builtins.open", mocker.mock_open(read_data="foo\t1\nbar\t5\n"))
    result = calculate_genome_size("path")
    assert result == 6


def test_calculate_average_coverage():
    result = calculate_average_coverage(
        genome_size=10, aligned_read_count=3, read_length=3
    )
    assert isinstance(result, dict)
    assert result["average_coverage"] == pytest.approx(0.9)


def test_make_qc_record():
    metric_1 = ("foo", {"bar": "baz"})
    metric_2 = ("qux", {"quux": "corge"})
    result = make_qc_record([metric_1, metric_2])
    assert result.to_ordered_dict["foo"] == {"bar": "baz"}
    assert result["qux"] == {"quux": "corge"}
    assert isinstance(result, QCMetricRecord)
