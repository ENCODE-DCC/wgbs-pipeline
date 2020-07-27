"""
Calculates average coverage from a bam file. The formula for this is given by this:
Coverage = (Aligned read counts * read length ) / total genome size

The total genome size is obtained by summing the values in the chromosome sizes. The
aligned read counts is the number of reads in the bam. The read length is obtained from
samtools-stats `RL`, see http://www.htslib.org/doc/samtools-stats.html for details.
"""

import argparse
import tempfile
from typing import Any, Dict, List, Tuple

import pysam
from qc_utils import QCMetric, QCMetricRecord
from qc_utils.parsers import parse_samtools_stats


def main() -> None:  # pragma: no cover
    """
    We already integration test this in the WDL tests, no need for Python test coverage.
    """
    parser = get_parser()
    args = parser.parse_args()
    samtools_stats = get_samtools_stats(args.bam, args.threads)
    genome_size = calculate_genome_size(args.chrom_sizes)
    average_coverage = calculate_average_coverage(
        genome_size=genome_size,
        aligned_read_count=samtools_stats["reads mapped"],
        read_length=samtools_stats["average length"],
    )
    qc_record = make_qc_record(
        [("samtools_stats", samtools_stats), ("average_coverage", average_coverage)]
    )
    qc_record.save(args.outfile)


def get_samtools_stats(bam_path: str, threads: int) -> Dict[str, Any]:
    """
    The buffer must be flushed first before it is read again by the parser, see
    https://stackoverflow.com/questions/46004774/python-namedtemporaryfile-appears-empty-even-after-data-is-written
    """
    samtools_stats_stdout = pysam.stats("--threads", str(threads), bam_path)
    with tempfile.NamedTemporaryFile(mode="w") as samtools_stats_file:
        samtools_stats_file.write(samtools_stats_stdout)
        samtools_stats_file.flush()
        samtools_stats = parse_samtools_stats(samtools_stats_file.name)
    return samtools_stats


def calculate_genome_size(chrom_sizes_path: str) -> int:
    size = 0
    with open(chrom_sizes_path) as f:
        for line in f:
            size += int(line.split()[1])
    return size


def calculate_average_coverage(
    genome_size: int, aligned_read_count: int, read_length: int
) -> Dict[str, float]:
    average_coverage = (aligned_read_count * read_length) / genome_size
    return {"average_coverage": average_coverage}


def make_qc_record(qcs: List[Tuple[str, Dict[str, Any]]]) -> QCMetricRecord:
    qc_record = QCMetricRecord()
    for name, qc in qcs:
        metric = QCMetric(name, qc)
        qc_record.add(metric)
    return qc_record


def get_parser() -> argparse.ArgumentParser:  # pragma: nocover
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", help="Bam file to compute coverage on", required=True)
    parser.add_argument("--chrom-sizes", help="Chromosome sizes file", required=True)
    parser.add_argument(
        "--threads", help="Number of threads for Samtools", required=True
    )
    parser.add_argument("--outfile", required=True)
    return parser


if __name__ == "__main__":
    main()
