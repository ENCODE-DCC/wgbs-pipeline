import argparse
import csv
from pathlib import Path
from typing import List

CSV_FIELDS = ["Barcode", "Name", "Dataset", "File1", "File2"]


def main():
    parser = get_parser()
    args = parser.parse_args()
    rows = process(args)
    with open("{}_metadata.csv".format(args.sample_names[0]), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(rows)


def process(args) -> List[List[str]]:
    """
    If only one file is specified (single-ended) then the File2 column is omitted.
    """
    files = read_file_tsv(args.files)
    if len(args.sample_names) != len(files):
        raise ValueError("Number of samples must match number of file groups")
    num_files = len(files[0])
    if num_files == 1:
        header = CSV_FIELDS[:-1]
    else:
        header = CSV_FIELDS
    to_write = [header]
    for sample_name, sample_files in zip(args.sample_names, files):
        output = ["sample_{}".format(sample_name), sample_name, sample_name]
        output.extend(sample_files)
        to_write.append(output)
    return to_write


def read_file_tsv(path: str) -> List[List[str]]:
    """
    Read a tsv where each row should contain one or two fastq files, and preprocess them
    by truncating paths to the basename.
    """
    output = []
    with open(path, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) not in (1, 2):
                raise ValueError("Expected at most two fastq files as input.")
            output.append([Path(i).name for i in row])
    return output


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-n", "--sample-names", nargs="+", help="A list of sample names", required=True
    )
    parser.add_argument(
        "--files",
        help="A tsv where each row is a fastq file or a pair of fastqs",
        required=True,
    )
    parser.add_argument(
        "-b",
        "--barcode-prefix",
        help="Prefix to prepend to sample names when making barcodes",
        default="sample_",
    )
    return parser


if __name__ == "__main__":
    main()
