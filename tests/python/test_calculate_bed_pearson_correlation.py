from contextlib import suppress as does_not_raise
from io import BytesIO

import pandas as pd
import pytest

from wgbs_pipeline.calculate_bed_pearson_correlation import (
    DF_COLUMN_NAMES,
    calculate_pearson,
    get_parser,
    load_bedmethyl,
    make_pearson_qc,
)


def test_load_bedmethyl():
    data = BytesIO(
        (
            b'track name="ENCSR156JXJ" description="ENCSR156JXJ" visibility=2 itemRgb="'
            b'On"\nchr1\t10468\t10469\t"ENCSR156JXJ"\t1\t+\t10468\t10469\t255,0,0\t1\t1'
            b"00\tCG\tCG\t2"
        )
    )
    df = load_bedmethyl(data)
    assert df.shape[0] == 1
    expected = ("chr1", 10468, 10489, 1, 100)
    for key, expected in zip(DF_COLUMN_NAMES, expected):
        assert df.at[0, key] == expected


def test_calculate_pearson():
    df1 = pd.DataFrame.from_records(
        [
            ("chr1", "1", "2", "14", "57"),
            ("chr1", "3", "4", "12", "100"),
            ("chr1", "5", "6", "3", "100"),
        ],
        columns=DF_COLUMN_NAMES,
    )
    df2 = pd.DataFrame.from_records(
        [
            ("chr1", "1", "2", "16", "87"),
            ("chr1", "3", "4", "12", "33"),
            ("chr1", "5", "6", "3", "100"),
            ("chr1", "7", "8", "3", "100"),
        ],
        columns=DF_COLUMN_NAMES,
    )
    pearson = calculate_pearson(df1, df2)
    assert pearson == 0.12


def test_make_pearson_qc():
    qc = make_pearson_qc(0.33)
    assert dict(qc.to_ordered_dict()) == {
        "pearson_correlation": {"pearson_correlation": 0.33}
    }


@pytest.mark.parametrize(
    "condition,args",
    [
        (does_not_raise(), ["--bedmethyls", "foo", "bar", "--outfile", "baz"]),
        (pytest.raises(SystemExit), ["--bedmethyls", "foo", "--outfile", "baz"]),
    ],
)
def test_get_parser(condition, args):
    parser = get_parser()
    with condition:
        parser.parse_args(args)
