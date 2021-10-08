import argparse
import csv
import json
from pathlib import Path
from typing import List

CSV_FIELDS = ["Barcode", "Name", "Dataset", "File1", "File2"]


def main():
    parser = get_parser()
    args = parser.parse_args()
    rows = process(args)
    with open(args.outfile, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(rows)


def process(args) -> List[List[str]]:
    """
    If only one file is specified (single-ended) then the File1 and File2 columns are
    supplanted by a single File column.
    """
    files = read_file_json(args.files)
    if len(args.sample_names) != len(files):
        raise ValueError("Number of samples must match number of file groups")
    num_files = len(files[0][0])
    if num_files == 1:
        header = CSV_FIELDS[:-2] + ["File"]
    else:
        header = CSV_FIELDS
    to_write = [header]
    technical_replicate_number = 0
    for sample_name, biological_replicate in zip(args.sample_names, files):
        constant_output = ["{}{}".format(args.barcode_prefix, sample_name), sample_name]
        for technical_replicate in biological_replicate:
            output = constant_output + [str(technical_replicate_number)] + technical_replicate
            technical_replicate_number += 1
            to_write.append(output)
    return to_write


def read_file_json(path: str) -> List[List[List[str]]]:
    """
    Read a JSON file of thrice nested arrays where all of the innermost arrays should
    contain the same number of fastq files (1 or 2), and preprocess them by truncating
    paths to the basename. We preserve the original JSON structure to avoid collapsing
    the biological replicate information.
    """
    output = []
    with open(path) as f:
        data = json.load(f)
        num_files = len(data[0][0])
        if num_files > 2:
            raise ValueError(
                "At most two files may be specified per techincal replicate, found "
                f"{num_files}"
            )
        for biological_replicate in data:
            paths_by_tech_rep = []
            for technical_replicate in biological_replicate:
                current_length = len(technical_replicate)
                if current_length != num_files:
                    raise ValueError(
                        (
                            "Found mixed paired and single ended libraries, the first "
                            f"library has {num_files} files but the current library "
                            f"contains {current_length}"
                        )
                    )
                paths_by_tech_rep.append([Path(i).name for i in technical_replicate])
            output.append(paths_by_tech_rep)
    return output


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-n", "--sample-names", nargs="+", help="A list of sample names", required=True
    )
    parser.add_argument(
        "-f",
        "--files",
        help="A JSON file containing an array of arrays of arrays",
        required=True,
    )
    parser.add_argument(
        "-b",
        "--barcode-prefix",
        help="Prefix to prepend to sample names when making barcodes",
        default="sample_",
    ),
    parser.add_argument(
        "-o",
        "--outfile",
        help="Name of csv file to write generated metadata to",
        default="metadata.csv",
    )
    return parser


if __name__ == "__main__":
    main()
