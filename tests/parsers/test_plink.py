import pytest

from cgr_gwas_qc.parsers import plink


@pytest.mark.real_data
def test_read_het(real_data_cache):
    filename = real_data_cache / "legacy_outputs/plink_start/samples_start.het"
    df = plink.read_het(filename)
    assert "ID" == df.index.name
    assert all(["O_HOM", "E_HOM", "N_NM", "F"] == df.columns.values)
