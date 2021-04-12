import pytest

from cgr_gwas_qc.cli.config import get_number_snps, get_output_pattern
from cgr_gwas_qc.testing.data import FakeData


def test_get_number_snps():
    bpm_file = FakeData._data_path / FakeData._illumina_manifest_file
    assert 2000 == get_number_snps(bpm_file)


@pytest.mark.parametrize(
    "sample_sheet_stem,output_pattern",
    [
        (
            "SR0001-001_1_AnalysisManifest_0000000000000",
            "{prefix}/SR0001-001_1_{file_type}_0000000000000.{ext}",
        ),
        ("test_AnalysisManifest_0123", "{prefix}/test_{file_type}_0123.{ext}"),
        ("test_0123", "{prefix}/{file_type}.{ext}"),
    ],
)
def test_get_output_pattern(sample_sheet_stem, output_pattern):
    assert output_pattern == get_output_pattern(sample_sheet_stem)
