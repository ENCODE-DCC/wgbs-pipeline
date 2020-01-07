import argparse
import csv
from pathlib import Path
from typing import List

CSV_FIELDS = ["Barcode", "Name", "Dataset", "File1", "File2"]


def main():
    parser = get_parser()
    args = parser.parse_args()
    rows = process(args)
    with open("{}_metadata.csv".format(args.sample_name), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(rows)


def process(args) -> List[List[str]]:
    """
    If only one file is specified (single-ended) then the File2 column is omitted.
    """
    num_files = len(args.files)
    if num_files > 2:
        raise ValueError("Expected at most two fastq files as input.")
    output = ["sample_{}".format(args.sample_name), args.sample_name, args.sample_name]
    file_basenames = [Path(i).name for i in args.files]
    if num_files == 1:
        header = CSV_FIELDS[:-1]
    else:
        header = CSV_FIELDS
    output.extend(file_basenames)
    return [header, output]


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--sample-name", help="sample name", required=True)
    parser.add_argument(
        "--files",
        nargs="+",
        help="paths to fastq files, must be at least one and at most two",
        required=True,
    )
    return parser


if __name__ == "__main__":
    main()
