import pytest

from cgr_gwas_qc.cli.config import get_number_snps, get_output_pattern


def test_get_number_snps(bpm_file):
    assert 2000 == get_number_snps(bpm_file)


@pytest.mark.parametrize(
    "sample_sheet_name,output_pattern",
    [
        ("test_AnalysisManifest_0123.csv", "{prefix}/test_{file_type}_0123.{ext}"),
        ("test_0123.csv", "{prefix}/{file_type}.{ext}"),
    ],
)
def test_get_output_pattern(sample_sheet_name, output_pattern):
    assert output_pattern == get_output_pattern(sample_sheet_name)
