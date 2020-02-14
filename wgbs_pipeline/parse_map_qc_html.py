import argparse
import json
import logging
from typing import Dict, List, Union

from bs4 import BeautifulSoup


def main():
    parser = get_parser()
    args = parser.parse_args()
    with open(args.infile) as f:
        html_doc = f.read()
    qc = parse_map_qc_html(html_doc)
    with open(args.outfile, "w") as f:
        json.dump(qc, f)


def parse_map_qc_html(html_doc: str) -> Dict[str, Union[float, int]]:
    qc: Dict[str, Union[float, int]] = {}
    soup = BeautifulSoup(html_doc, "html.parser")

    map_stats_table = find_table(soup, name="Mapping Stats (Reads)")
    for row in map_stats_table:
        qc_key = qc_key_from_tag(row.td)
        qc_value: Union[int, float] = int_from_tag(row.contents[3])
        qc_pct_key = "pct_" + qc_key
        qc_pct_value = percent_from_tag(row.contents[5])
        qc[qc_key] = qc_value
        qc[qc_pct_key] = qc_pct_value

    fragment_uniqueness_table = find_table(soup, name="Uniqueness (Fragments)")
    for row in fragment_uniqueness_table:
        qc_key = qc_key_from_tag(row.td)
        if qc_key == "average_unique":
            qc_key = "pct_unique_fragments"
            qc_value = percent_from_tag(row.contents[3])
        else:
            qc_value = int_from_tag(row.contents[3])
        qc[qc_key] = qc_value

    bisulfite_conversion_rate_table = find_table(soup, name="Bisulfite Conversion Rate")
    for row in bisulfite_conversion_rate_table:
        qc_key = qc_key_from_tag(row.td)
        try:
            qc_value = float(string_from_tag(row.contents[3]))
        except ValueError:
            logging.warning(
                f"Could not parse qc value for {qc_key} into float", exc_info=True
            )
            continue
        qc[qc_key] = qc_value

    correct_pairs_table = find_table(soup, name="Correct Pairs")
    for row in correct_pairs_table:
        qc_key = qc_key_from_tag(row.td)
        qc_value = int_from_tag(row.contents[3])
        qc[qc_key] = qc_value

    return qc


def find_table(soup: BeautifulSoup, name: str) -> List[BeautifulSoup]:
    """
    Extract the table in the HTML by name and return its child elements, ignoring
    nonsensical rows like "\n" and " ". Also ignores <th> elements (table headers) since
    they aren't useful.
    """
    output = []
    table = soup.find(
        lambda tag: tag.name == "h1" and tag.string == f" {name} "
    ).next_sibling.next_sibling.contents
    for row in table:
        if row in ("\n", " ") or row.th is not None:
            continue
        output.append(row)
    return output


def percent_from_tag(percent: BeautifulSoup) -> float:
    """
    Takes a string like " 99.0 % " and converts it to a float in the range [0, 1]
    """
    return float(string_from_tag(percent).rstrip("%").strip()) / 100


def string_from_tag(tag: BeautifulSoup) -> str:
    return tag.string.strip()


def int_from_tag(tag: BeautifulSoup) -> int:
    return int(string_from_tag(tag))


def qc_key_from_tag(tag: BeautifulSoup) -> str:
    """
    Convert the row names in the HTML table to portal-style property names by converting
    to snake case. If prepend_pct is true, then return a tuple of names instead, with
    the second having "pct_" prepended to it.
    """
    return string_from_tag(tag).lower().replace(" ", "_")


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="Input HTML file", required=True)
    parser.add_argument(
        "-o", "--outfile", help="Name of output JSON file", required=True
    )
    return parser


if __name__ == "__main__":
    main()
