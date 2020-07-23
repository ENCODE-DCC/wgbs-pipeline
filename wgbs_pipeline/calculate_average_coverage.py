"""
Calculates average coverage from a bam file. The formula for this is given by this:
Coverage = (Aligned read counts * read length ) / total genome size

The total genome size is obtained by summing the values in the chromosome sizes. The
aligned read counts is the number of reads in the bam. The read length is obtained from
samtools-stats `RL`, see http://www.htslib.org/doc/samtools-stats.html for details.
"""

import argparse
from io import StringIO
from typing import Dict

import pysam
from qc_utils import QCMetric, QCMetricRecord
from qc_utils.parsers import parse_samtools_stats


def main() -> None:  # pragma: no cover
    """
    We already integration test this in the WDL tests, no need for Python test coverage.
    """
    parser = get_parser()
    args = parser.parse_args()
    samtools_stats_stdout = pysam.stats("--threads", args.threads, args.bam)
    samtools_stats_buffer = StringIO(samtools_stats_stdout)
    samtools_stats = parse_samtools_stats(samtools_stats_buffer)
    average_coverage = calculate_average_coverage(1, 2, 3)
    qc_record = QCMetricRecord()
    samtools_stats_metric = QCMetric("samtools_stats", samtools_stats)
    average_coverage_metric = QCMetric("average_coverage", average_coverage)
    qc_record.add(samtools_stats_metric)
    qc_record.add(average_coverage_metric)
    qc_record.save(args.outfile)


def calculate_genome_size(chrom_sizes_path: str):
    pass


def calculate_read_length(bam_path: str):
    pass


def calculate_average_coverage(
    genome_size: int, aligned_read_counts: int, read_length: int
) -> Dict[str, float]:
    return {"average_coverage": 1.0}


def get_parser() -> argparse.ArgumentParser:
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
