import pandas as pd
import pytest
from pandas.testing import assert_series_equal

from cgr_gwas_qc.testing.data import RealData
from cgr_gwas_qc.workflow.scripts.snp_qc_table import build_table


@pytest.fixture
def thousand_genomes_fake(tmp_path):
    csv = tmp_path / "test.csv"

    csv.write_text(
        "\n".join(
            [
                "array_id,thousand_genomes_id",
                "2010-08-Y-1110,rs12345Fake",
                "2010-08-Y-1111,rs12346Fake",
                "2010-08-Y-1221,rs12376Fake",
                "MitoA16163G,rs99999Fake",
            ]
        )
    )
    return csv


@pytest.mark.real_data
def test_snp_qc_table(thousand_genomes_fake):
    data_cache = RealData()
    initial = data_cache / "production_outputs/plink_start/samples_start.lmiss"
    cr1 = data_cache / "production_outputs/plink_filter_call_rate_1/samples_filter1.lmiss"
    cr2 = data_cache / "production_outputs/plink_filter_call_rate_2/samples_filter2.lmiss"
    df = build_table(initial, cr1, cr2, thousand_genomes_fake)

    assert_series_equal(
        pd.read_csv(thousand_genomes_fake).set_index("array_id").squeeze(),
        df[["array_snp_id", "thousand_genomes_snp_id"]]
        .dropna()
        .set_index("array_snp_id")
        .squeeze(),
        check_names=False,
    )
    assert_series_equal(df.snp_cr1.isnull(), df.snp_cr1_removed, check_names=False)
    assert_series_equal(df.snp_cr2.isnull(), df.snp_cr2_removed, check_names=False)
