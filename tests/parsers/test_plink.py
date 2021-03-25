import pytest

from cgr_gwas_qc.parsers import plink
from cgr_gwas_qc.testing.data import RealData


@pytest.mark.real_data
def test_read_het():
    filename = RealData() / "production_outputs/plink_start/samples_start.het"
    df = plink.read_het(filename)
    assert "ID" == df.index.name
    assert all(["O_HOM", "E_HOM", "N_NM", "F"] == df.columns.values)
