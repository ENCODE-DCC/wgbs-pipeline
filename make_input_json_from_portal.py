import argparse
import json
from pathlib import Path
from urllib.parse import urljoin

try:
    import requests
except ImportError:
    print("Make sure to install requests")
    raise

PORTAL_URL = "https://www.encodeproject.org"

REFERENCE_FILES = {
    "GRCh38": {
        "reference": "/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz",
        "extra_reference": "/files/lambda.fa/@@download/lambda.fa.fasta.gz",
        "indexed_reference": urljoin(
            PORTAL_URL, "/files/ENCFF603ORM/@@download/ENCFF603ORM.tar.gz"
        ),
        "indexed_contig_sizes": "/files/ENCFF792NJK/@@download/ENCFF792NJK.tsv",
    },
    "mm10": {
        "reference": "/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz",
        "extra_reference": "/files/lambda.fa/@@download/lambda.fa.fasta.gz",
        "indexed_reference": urljoin(
            PORTAL_URL, "/files/ENCFF708YIO/@@download/ENCFF708YIO.tar.gz"
        ),
        "indexed_contig_sizes": "/files/mm10_no_alt.chrom.sizes/@@download/mm10_no_alt.chrom.sizes.tsv",
    },
}
ALLOWED_STATUSES = ("released", "in progress")


def main():
    parser = get_parser()
    args = parser.parse_args()
    auth = _read_auth_from_file(args.keypair_file)
    experiment = get_experiment(args.accession, auth=auth)
    fastqs = _get_fastqs_from_experiment(experiment)
    assembly_name = _get_assembly_from_experiment(experiment)
    input_json = _get_input_json(fastqs=fastqs, assembly_name=assembly_name)
    outfile = args.outfile or "{}.json".format(args.accession)
    _write_json_to_file(input_json, outfile)


def get_experiment(accession, auth=None):
    response = requests.get(
        urljoin(PORTAL_URL, accession),
        auth=auth,
        headers={"Accept": "application/json"},
    )
    response.raise_for_status()
    return response.json()


def _get_fastqs_from_experiment(experiment):
    fastq_pairs_by_replicate = {}
    for file in experiment["files"]:
        if file["file_format"] == "fastq" and file["status"] in ALLOWED_STATUSES:
            biological_replicate = file["biological_replicates"][0]
            paired_with = file.get("paired_with")
            if paired_with is not None:
                if file["paired_end"] == "2":
                    continue
                paired_with_file = [
                    f for f in experiment["files"] if f["@id"] == paired_with
                ][0]
                fastq_pair = [
                    urljoin(PORTAL_URL, file["href"]),
                    urljoin(PORTAL_URL, paired_with_file["href"]),
                ]
            else:
                fastq_pair = [urljoin(PORTAL_URL, file["href"])]
            replicate_fastqs = fastq_pairs_by_replicate.get(biological_replicate)
            if replicate_fastqs is None:
                fastq_pairs_by_replicate[biological_replicate] = []
            fastq_pairs_by_replicate[biological_replicate].append(fastq_pair)
    output = [replicate for replicate in fastq_pairs_by_replicate.values()]
    return output


def _get_input_json(fastqs, assembly_name):
    input_json = {
        "wgbs.fastqs": fastqs,
        "wgbs.reference": REFERENCE_FILES[assembly_name]["reference"],
        "wgbs.extra_reference": REFERENCE_FILES[assembly_name]["extra_reference"],
        "wgbs.indexed_reference": REFERENCE_FILES[assembly_name]["indexed_reference"],
        "wgbs.indexed_contig_sizes": REFERENCE_FILES[assembly_name][
            "indexed_contig_sizes"
        ],
    }
    return input_json


def _write_json_to_file(data, outfile):
    Path(outfile).write_text(json.dumps(data, indent=2, sort_keys=True))


def _read_auth_from_file(keypair_file):
    keypair_path = Path(keypair_file).expanduser()
    if keypair_path.exists():
        data = json.loads(keypair_path.read_text())
        return (data["submit"]["key"], data["submit"]["secret"])
    else:
        return None


def _get_assembly_from_experiment(experiment):
    organism_name = experiment["replicates"][0]["library"]["biosample"]["organism"][
        "name"
    ]
    if organism_name == "human":
        return "GRCh38"
    if organism_name == "mouse":
        return "mm10"
    raise ValueError(f"Could not determine assembly for organism {organism_name}")


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "accession", help="Accession of portal experiment to generate input for"
    )
    parser.add_argument(
        "-o",
        "--outfile",
        help="Name of file to write output JSON to, defaults to [ACCESSION].json",
    )
    parser.add_argument(
        "--keypair-file", help="Path to keypairs.json", default="~/keypairs.json"
    )
    return parser


if __name__ == "__main__":
    main()
