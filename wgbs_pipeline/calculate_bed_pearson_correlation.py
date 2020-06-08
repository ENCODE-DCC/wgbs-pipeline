"""
Calculates Pearson correlation of sites with ≥10X coverage between CpG bedmethyl files.
Inspired by https://github.com/ENCODE-DCC/mirna-seq-pipeline/blob/dev/src/calculate_correlation.py
"""

import argparse
import json

import pandas as pd
from qc_utils import QCMetric, QCMetricRecord

DF_COLUMN_NAMES = ("chrom", "start", "end", "coverage", "methylation_percentage")


def main() -> None:
    parser = get_parser()
    args = parser.parse_args()
    bedmethyl1 = load_bedmethyl(args.bedmethyls[0])
    bedmethyl2 = load_bedmethyl(args.bedmethyls[1])
    pearson_correlation = calculate_pearson(bedmethyl1, bedmethyl2)
    qc_record = make_pearson_qc(pearson_correlation)

    # Should be able to do something like qc_record.save(args.outfile) in the future.
    with open(args.outfile, "w") as fp:
        json.dump(qc_record.to_ordered_dict(), fp)


def load_bedmethyl(path: str) -> pd.DataFrame:
    """
    Loads position information (columns 0-2), coverage (column 10), and methylation
    percentage of the input bedMethyl file, for file specifications see
    https://www.encodeproject.org/data-standards/wgbs/ .
    """
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        skiprows=1,
        usecols=[0, 1, 2, 9, 10],
        names=DF_COLUMN_NAMES,
    )
    return df


def calculate_pearson(bedmethyl1: pd.DataFrame, bedmethyl2: pd.DataFrame) -> float:
    """
    Filters both bedfiles for entries where coverage is GTE 10, then does a join to get
    the positions to line up, and calculated the correlation.
    """
    bedmethyl1 = bedmethyl1.loc([bedmethyl1["coverage"] >= 10])
    bedmethyl2 = bedmethyl2.loc([bedmethyl2["coverage"] >= 10])
    pearson_correlation = bedmethyl1.corr(bedmethyl2, method="pearson")
    return pearson_correlation


def make_pearson_qc(pearson_correlation: float) -> QCMetricRecord:
    qc_record = QCMetricRecord()
    pearson_metric = QCMetric(
        "pearson_correlation", {"pearson_correlation": pearson_correlation}
    )
    qc_record.add(pearson_metric)
    return qc_record


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bedmethyls", nargs=2, help="CpG bedMethyl files to correlate.", required=True
    )
    parser.add_argument("-o", "--outfile", required=True)
    return parser


if __name__ == "__main__":
    main()
