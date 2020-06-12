import pytest


@pytest.mark.workflow("test_wgbs_two_reps")
def test_wgbs_two_reps_check_pearson_qc(workflow_dir, test_data_dir):
    with open(test_data_dir / "bed_pearson_correlation_qc.json") as f:
        expected = f.read()

    result_path = next(
        (workflow_dir / "test-output").glob(
            "wgbs/*/call-calculate_bed_pearson_correlation/execution/bed_pearson_correlation_qc.json"
        )
    )
    with open(result_path) as f:
        result = f.read()
    assert result == expected
