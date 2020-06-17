import json
import math

import pytest


@pytest.mark.workflow("test_wgbs_two_reps")
def test_wgbs_two_reps_check_pearson_qc(workflow_dir, test_data_dir):
    with open(test_data_dir / "bed_pearson_correlation_qc.json") as f:
        expected_qc = json.load(f)

    result_path = next(
        (workflow_dir / "test-output").glob(
            "wgbs/*/call-calculate_bed_pearson_correlation/execution/bed_pearson_correlation_qc.json"
        )
    )
    with open(result_path) as f:
        result_qc = json.load(f)

    result = result_qc["pearson_correlation"]["pearson_correlation"]
    expected = expected_qc["pearson_correlation"]["pearson_correlation"]
    assert math.isclose(result, expected, rel_tol=1e-5)
