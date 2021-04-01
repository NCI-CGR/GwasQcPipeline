from io import StringIO
from textwrap import dedent

import pandas as pd
import pytest

from cgr_gwas_qc.models.config import Config
from cgr_gwas_qc.testing.data import RealData
from cgr_gwas_qc.workflow.scripts.population_qc_table import _expand_related, main


def test_expand_related():
    relatives = dedent(
        """
        QC_Family_ID,relatives
        fam1,S001|S002
        """
    )

    df = pd.concat(_expand_related(StringIO(relatives), "test"), ignore_index=True)

    assert (2, 3) == df.shape
    assert ["S001", "S002"] == sorted(df.Subject_ID)
    assert {"S001|S002"} == set(df.relatives)


@pytest.mark.real_data
def test_main(tmp_path, real_config: Config):
    pca = RealData() / "production_outputs/pca/EUR_subjects.eigenvec"
    het = RealData() / "production_outputs/autosomal_heterozygosity/EUR_subjects_qc.het"
    relatives = dedent(
        """
        QC_Family_ID,relatives
        fam1,W-118174-2|D-144974-4
        fam2,F-094558-7|F-141988-8
        """
    )

    main(
        StringIO(relatives),
        pca,
        het,
        "test",
        real_config.software_params.autosomal_het_threshold,
        tmp_path / "test.csv",
    )

    df = pd.read_csv(tmp_path / "test.csv")

    assert df.query("QC_Family_ID == 'test_fam1'").shape[0] == 2
    assert df.query("QC_Family_ID == 'test_fam2'").shape[0] == 2
    assert df.dropna(how="any").shape[0] == 4
