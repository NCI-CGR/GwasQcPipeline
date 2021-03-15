import numpy as np
import pytest

from cgr_gwas_qc.parsers.eigensoft import Eigenvec
from cgr_gwas_qc.testing.data import RealData


@pytest.fixture(scope="module")
def eigenvec():
    filename = RealData() / "production_outputs/pca/EUR_subjects.eigenvec"
    return Eigenvec(filename)


def test_Eigenvec_componenets(eigenvec: Eigenvec):
    cols = ["PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"]
    assert "ID" == eigenvec.components.index.name
    assert cols == eigenvec.components.columns.tolist()
    assert all(eigenvec.components.dtypes[1:] == np.float64)  # all PCs should be floats


def test_Eigenvec_values(eigenvec: Eigenvec):
    assert eigenvec.values.shape[0] == 10
    assert eigenvec.values.dtype == np.float64
