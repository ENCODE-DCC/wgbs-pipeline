import argparse
from pathlib import Path
from typing import List


def main():
    parser = get_parser()
    args = parser.parse_args()
    conf = make_conf(args)
    write_conf(args.outfile, conf)


def make_conf(args: argparse.Namespace) -> List[str]:
    """
    The reference and extra reference files are localized in the reference directory at
    the WDL level. Therefore, we only need the basename. However, for the include file,
    the actual path must be provided in order for gemBS to recognize it.
    """
    conf = [
        f"reference = reference/{Path(args.reference).name}",
        f"extra_references = reference/{Path(args.extra_reference).name}",
        "index_dir = indexes",
        "base = .",
        "sequence_dir = ${base}/fastq/@SAMPLE",
        "bam_dir = ${base}/mapping/@BARCODE",
        "bcf_dir = ${base}/calls/@BARCODE",
        "extract_dir = ${base}/extract/@BARCODE",
        "report_dir = ${base}/report",
        f"threads = {args.num_threads}",
        f"jobs = {args.num_jobs}",
    ]

    if args.underconversion_sequence or args.include_file:
        conf.append("[mapping]")
    if args.underconversion_sequence:
        conf.append(f"underconversion_sequence = {args.underconversion_sequence}")
    if args.include_file:
        conf.append(f"include {args.include_file}")

    return conf


def write_conf(outfile: str, conf: List[str]):
    with open(outfile, "w") as f:
        for line in conf:
            f.write(line + "\n")


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-t", "--num-threads", help="Number of threads for gemBS", type=int, default=8
    )
    parser.add_argument(
        "-j", "--num-jobs", help="Number of jobs for gemBS", type=int, default=3
    )
    parser.add_argument(
        "-e",
        "--extra-reference",
        help="File name of extra reference for control sequence",
        required=True,
    )
    parser.add_argument(
        "-r", "--reference", help="File name of reference fasta", required=True
    )
    parser.add_argument(
        "-u",
        "--underconversion-sequence",
        help="Name of contig in extra reference to use as underconversion control",
    )
    parser.add_argument(
        "-i", "--include-file", help="Name of additional conf file to include in "
    )
    parser.add_argument(
        "-o",
        "--outfile",
        help="Name of file to write generated conf to",
        default="gembs.conf",
    )
    return parser


if __name__ == "__main__":
    main()
