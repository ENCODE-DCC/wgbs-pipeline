import builtins
from contextlib import suppress as does_not_raise

import pytest
from bs4 import BeautifulSoup

from wgbs_pipeline.parse_map_qc_html import (
    find_table,
    get_parser,
    int_from_tag,
    main,
    parse_map_qc_html,
    percent_from_tag,
    qc_key_from_tag,
    string_from_tag,
)


def test_main(mocker):
    """
    The assert looks a little wonky here. The first index extracts the args of the call,
    and the second extracts the first positional arg.
    """
    mocker.patch("builtins.open", mocker.mock_open(read_data=HTML_DOC))
    testargs = ["prog", "-i", "infile.html", "-o", "outfile.json"]
    mocker.patch("sys.argv", testargs)
    main()
    assert builtins.open.call_args[0][0] == "outfile.json"


def test_parse_map_qc_html():
    """
    Need to use pytest.approx here due to floating point comparisons
    """
    qc = parse_map_qc_html(HTML_DOC)
    assert qc == pytest.approx(
        {
            "sequenced_reads": 400000,
            "pct_sequenced_reads": 1.0,
            "general_reads": 358466,
            "pct_general_reads": 0.8962,
            "reads_in_control_sequences": 0,
            "pct_reads_in_control_sequences": 0.0,
            "reads_under_conversion_control": 12349,
            "pct_reads_under_conversion_control": 0.0309,
            "reads_over_conversion_control": 0,
            "pct_reads_over_conversion_control": 0.0,
            "unmapped_reads": 29185,
            "pct_unmapped_reads": 0.073,
            "bisulfite_reads_c2t": 192236,
            "pct_bisulfite_reads_c2t": 0.4806,
            "bisulfite_reads_g2a": 178579,
            "pct_bisulfite_reads_g2a": 0.4464,
            "unique_fragments": 156841,
            "pct_unique_fragments": 0.7842,
            "conversion_rate": 0.9986559741624927,
            "correct_pairs": 179763,
        }
    )


def test_parse_map_qc_html_single_ended(single_ended_map_qc_html):
    """
    Single ended data is missing pair information. Need to use pytest.approx here due to
    floating point comparisons
    """
    qc = parse_map_qc_html(single_ended_map_qc_html)
    assert qc == pytest.approx(
        {
            "bisulfite_reads_c2t": 427911721,
            "bisulfite_reads_g2a": 410279399,
            "conversion_rate": 0.993987188263205,
            "general_reads": 834427445,
            "pct_bisulfite_reads_c2t": 0.3791,
            "pct_bisulfite_reads_g2a": 0.3635,
            "pct_general_reads": 0.7392,
            "pct_reads_in_control_sequences": 0.0,
            "pct_reads_over_conversion_control": 0.0,
            "pct_reads_under_conversion_control": 0.0033,
            "pct_sequenced_reads": 1.0,
            "pct_unique_fragments": 0.7465999999999999,
            "pct_unmapped_reads": 0.25739999999999996,
            "reads_in_control_sequences": 0,
            "reads_over_conversion_control": 0,
            "reads_under_conversion_control": 3763675,
            "sequenced_reads": 1128751544,
            "unique_fragments": 625822824,
            "unmapped_reads": 290560424,
        }
    )


def test_string_from_tag():
    tag = BeautifulSoup("<title> Foo </title>", "html.parser")
    result = string_from_tag(tag)
    assert result == "Foo"


def test_percent_from_tag():
    tag = BeautifulSoup("<title> 50.0 % </title>", "html.parser")
    result = percent_from_tag(tag)
    assert result == 0.5


def test_int_from_tag():
    tag = BeautifulSoup("<title> 3 </title>", "html.parser")
    result = int_from_tag(tag)
    assert result == 3


def test_qc_key_from_tag():
    tag = BeautifulSoup("<title> Mapped Reads </title>", "html.parser")
    result = qc_key_from_tag(tag)
    assert result == "mapped_reads"


def test_find_table(html_table):
    soup = BeautifulSoup(html_table, "html.parser")
    table = find_table(soup, name="Mapping Stats (Reads)")
    assert table[0].td.string == " Sequenced Reads "
    assert "\n" not in table
    assert all(i.th is None for i in table)


def test_find_table_not_found_returns_none(html_table):
    soup = BeautifulSoup(html_table, "html.parser")
    table = find_table(soup, name="foo")
    assert table is None


@pytest.mark.parametrize(
    "args,condition",
    [
        (["-i", "infile", "-o", "outfile"], does_not_raise()),
        (["--infile", "foo.html"], pytest.raises(SystemExit)),
        (["--outfile", "out.json"], pytest.raises(SystemExit)),
    ],
)
def test_get_parser(args, condition):
    parser = get_parser()
    with condition:
        parser.parse_args(args)


HTML_DOC = """
<HTML>
 <HEAD>
 </HEAD>
 <BODY>
 <P id='path'> /ENCODE/sample_sample5 </P>
 <a class="link" href="ENCODE.html"><B>BACK</B></a> <br>
  <H1 id="title"> <U> SAMPLE sample_sample5 </U> </H1>
<BR><BR><BR>
<H1 id="section"> Mapping Stats (Reads) </H1>
  <TABLE id="hor-zebra">
   <TR>
    <TH scope="col">Concept</TH> <TH scope="col">Total Reads</TH> <TH scope="col">%</TH>
    <TH scope="col">Pair One Reads</TH> <TH scope="col">%</TH> <TH scope="col">Pair Two Reads</TH> <TH scope="col">%</TH>
   </TR>
   <TR class="odd">
   <TD> Sequenced Reads </TD> <TD> 400000 </TD> <TD> 100.00 % </TD> <TD> 200000 </TD> <TD> 100.00 % </TD> <TD> 200000 </TD> <TD> 100.00 % </TD>
   </TR>
   <TR>
   <TD> General Reads </TD> <TD> 358466 </TD> <TD> 89.62 % </TD> <TD> 179416 </TD> <TD> 89.71 % </TD> <TD> 179050 </TD> <TD> 89.53 % </TD>
   </TR>
   <TR class="odd">
   <TD> Reads in Control sequences </TD> <TD> 0 </TD> <TD> 0.00 % </TD> <TD> 0 </TD> <TD> 0.00 % </TD> <TD> 0 </TD> <TD> 0.00 % </TD>
   </TR>
   <TR>
   <TD> Reads under conversion control  </TD> <TD> 12349 </TD> <TD> 3.09 % </TD> <TD> 6179 </TD> <TD> 3.09 % </TD> <TD> 6170 </TD> <TD> 3.08 % </TD>
   </TR>
   <TR class="odd">
   <TD> Reads over conversion control </TD> <TD> 0 </TD> <TD> 0.00 % </TD> <TD> 0 </TD> <TD> 0.00 % </TD> <TD> 0 </TD> <TD> 0.00 % </TD>
   </TR>
   <TR>
   <TD> Unmapped reads </TD> <TD> 29185 </TD> <TD> 7.30 % </TD> <TD> 14405 </TD> <TD> 7.20 % </TD> <TD> 14780 </TD> <TD> 7.39 % </TD>
   </TR>
   <TR class="odd">
   <TD> Bisulfite_reads C2T </TD> <TD> 192236 </TD> <TD> 48.06 % </TD> <TD> 96211 </TD> <TD> 48.11 % </TD> <TD> 96025 </TD> <TD> 48.01 % </TD>
   </TR>
   <TR>
   <TD> Bisulfite_reads G2A </TD> <TD> 178579 </TD> <TD> 44.64 % </TD> <TD> 89384 </TD> <TD> 44.69 % </TD> <TD> 89195 </TD> <TD> 44.60 % </TD>
   </TR>
 </TABLE>
<BR><BR><BR>
<H1 id="section"> Uniqueness (Fragments) </H1>
  <TABLE id="green">
   <TR>
    <TH scope="col">Concept</TH> <TH scope="col">Value</TH>
   </TR>
  <TR>   <TD> Unique Fragments </TD> <TD> 156841 </TD>
  </TR>  <TR class="odd">   <TD> Average Unique </TD> <TD> 78.42 % </TD>
  </TR>
 </TABLE>
<BR><BR><BR>
<H1 id="section"> Bisulfite Conversion Rate </H1>
  <TABLE id="hor-zebra">
   <TR>
    <TH scope="col">Bisulfite Conversion Type</TH> <TH scope="col">Conversion Rate</TH>
   </TR>
  <TR>   <TD> Conversion Rate </TD> <TD> 0.9986559741624927 </TD>
  </TR>  <TR class="odd">   <TD> Over Conversion Rate </TD> <TD> NA </TD>
  </TR>
 </TABLE>
<BR><BR><BR>
<H1 id="section"> Correct Pairs </H1>
  <TABLE id="green">
   <TR>
    <TH scope="col">Concept</TH> <TH scope="col">Total Reads</TH>
   </TR>
  <TR>   <TD> Correct Pairs </TD> <TD> 179763 </TD>
   </TR>
 </TABLE>
 </BODY>
</HTML>
"""


@pytest.fixture
def html_table():
    table = """
    <HTML>
    <HEAD>
    </HEAD>
    <BODY>
    <H1 id="section"> Mapping Stats (Reads) </H1>
    <TABLE id="hor-zebra">
    <TR>
        <TH scope="col">Concept</TH> <TH scope="col">Total Reads</TH> <TH scope="col">%</TH>
        <TH scope="col">Pair One Reads</TH> <TH scope="col">%</TH> <TH scope="col">Pair Two Reads</TH> <TH scope="col">%</TH>
    </TR>
    <TR class="odd">
    <TD> Sequenced Reads </TD> <TD> 400000 </TD> <TD> 100.00 % </TD> <TD> 200000 </TD> <TD> 100.00 % </TD> <TD> 200000 </TD> <TD> 100.00 % </TD>
    </TR>
    </TABLE>
    </BODY>
    </HTML>
    """

    return table


@pytest.fixture
def single_ended_map_qc_html():
    html = """
    <HTML>
    <HEAD>
    <STYLE TYPE="text/css">
    <!--
    @import url("style.css");
    -->
    </STYLE>
    </HEAD>
    <BODY>

    <P id='path'> /ENCODE/sample_0 </P>

    <a class="link" href="ENCODE.html"><B>BACK</B></a> <br>
    <H1 id="title"> <U> SAMPLE sample_0 </U> </H1>
    <BR><BR><BR>
    <H1 id="section"> Mapping Stats (Reads) </H1>
    <TABLE id="hor-zebra">
    <TR>
        <TH scope="col">Concept</TH> <TH scope="col">Total Reads</TH> <TH scope="col">%</TH>
    </TR>
    <TR class="odd">
    <TD> Sequenced Reads </TD> <TD> 1128751544 </TD> <TD> 100.00 %</TD>
    </TR>
    <TR>
    <TD> General Reads </TD> <TD> 834427445 </TD> <TD> 73.92 %</TD>
    </TR>
    <TR class="odd">
    <TD> Reads in Control sequences </TD> <TD> 0 </TD> <TD> 0.00 %</TD>
    </TR>
    <TR>
    <TD> Reads under conversion control  </TD> <TD> 3763675 </TD> <TD> 0.33 %</TD>
    </TR>
    <TR class="odd">
    <TD> Reads over conversion control </TD> <TD> 0 </TD> <TD> 0.00 %</TD>
    </TR>
    <TR>
    <TD> Unmapped reads </TD> <TD> 290560424 </TD> <TD> 25.74 %</TD>
    </TR>
    <TR class="odd">
    <TD> Bisulfite_reads C2T </TD> <TD> 427911721 </TD> <TD> 37.91 %</TD>
    </TR>
    <TR>
    <TD> Bisulfite_reads G2A </TD> <TD> 410279399 </TD> <TD> 36.35 %</TD>
    </TR>
    </TABLE>
    <BR><BR><BR>
    <H1 id="section"> Uniqueness (Fragments) </H1>
    <TABLE id="green">
    <TR>
        <TH scope="col">Concept</TH> <TH scope="col">Value</TH>
    </TR>
    <TR>   <TD> Unique Fragments </TD> <TD> 625822824 </TD>
    </TR>  <TR class="odd">   <TD> Average Unique </TD> <TD> 74.66 % </TD>
    </TR>
    </TABLE>
    <BR><BR><BR>
    <H1 id="section"> Mapping Stats (Bases) </H1>
    <TABLE id="hor-zebra">
    <TR>
        <TH scope="col">Concept</TH> <TH scope="col">Total Bases</TH> <TH scope="col">%</TH>
    </TR>
    <TR class="odd">
    <TD> Base Counts Overall A </TD> <TD> 28551057960 </TD> <TD> 24.97 %</TD>
    </TR>
    <TR>
    <TD> Base Counts Overall C </TD> <TD> 1795515556 </TD> <TD> 1.57 %</TD>
    </TR>
    <TR class="odd">
    <TD> Base Counts Overall G </TD> <TD> 26644100723 </TD> <TD> 23.30 %</TD>
    </TR>
    <TR>
    <TD> Base Counts Overall T </TD> <TD> 44170064692 </TD> <TD> 38.63 %</TD>
    </TR>
    <TR class="odd">
    <TD> Base Counts Overall N </TD> <TD> 13191351907 </TD> <TD> 11.54 %</TD>
    </TR>
    </TABLE>
    <BR><BR><BR>
    <H1 id="section"> Bisulfite Conversion Rate </H1>
    <TABLE id="hor-zebra">
    <TR>
        <TH scope="col">Bisulfite Conversion Type</TH> <TH scope="col">Conversion Rate</TH>
    </TR>
    <TR>   <TD> Conversion Rate </TD> <TD> 0.993987188263205 </TD>
    </TR>  <TR class="odd">   <TD> Over Conversion Rate </TD> <TD> NA </TD>
    </TR>
    </TABLE>
    <BR><BR><BR>
    <H1 id="section"> Mapping Quality </H1>
    <TABLE id="hor-zebra">
    <TR class="odd"> <TH scope="col"> Mapping Quality Histogram</TH> </TR>
    <TR> <TD> <img src="sample_0.mapq.png" alt="sample_0.mapq.png"> </TD> </TR>
    </TABLE>
    <BR><BR><BR>
    <H1 id="section"> Mapping Lanes Reports </H1>
    <TABLE id="hor-zebra">
    <TR> <TH scope="col"> LANE REPORTS </TH> </TR>
    <TR class="odd"> <TD> <a class="link" href="0.html"> 0 </TD> </TR>
    </TABLE>
    </BODY>
    </HTML>
    """
    return html
