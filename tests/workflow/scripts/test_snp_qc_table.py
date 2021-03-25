from io import StringIO

import pandas as pd
import pytest
from pandas.testing import assert_series_equal

from cgr_gwas_qc.testing.data import RealData
from cgr_gwas_qc.workflow.scripts.snp_qc_table import add_call_rate_flags, aggregate_snp_data


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
def test_aggregate_snp_data(thousand_genomes_fake):
    # GIVEN: SNP Call Rate file and a fake mapping of thousand genome rsIDs
    data_cache = RealData()
    initial = data_cache / "production_outputs/plink_start/samples_start.lmiss"
    cr1 = data_cache / "production_outputs/plink_filter_call_rate_1/samples_filter1.lmiss"
    cr2 = data_cache / "production_outputs/plink_filter_call_rate_2/samples_filter2.lmiss"

    # WHEN: I aggregate and add call rate summaries
    df = aggregate_snp_data(initial, cr1, cr2, thousand_genomes_fake)

    # THEN: The thousand genomes file is correctly merged
    assert_series_equal(
        pd.read_csv(thousand_genomes_fake).set_index("array_id").squeeze(),
        df[["array_snp_id", "thousand_genomes_snp_id"]]
        .dropna()
        .set_index("array_snp_id")
        .squeeze(),
        check_names=False,
    )

    # Call rate columns should exist
    assert {"Call_Rate_Initial", "Call_Rate_1", "Call_Rate_2"}.issubset(set(df.columns))


def test_add_call_rate_flags():
    # GIVEN: A fake call rate table with each possible class
    df = pd.read_csv(
        StringIO(
            (
                "Call_Rate_Initial,Call_Rate_1,Call_Rate_2\n"
                ",,\n"  # Missing in initial
                "1,,\n"  # CR1 filtered
                "1,1,\n"  # CR2 filtered
                "1,1,1\n"  # Not filtered
            )
        )
    )

    # WHEN: I add summary flags
    add_call_rate_flags(df)

    # THEN: The number flagged should match the expectation
    assert df.is_cr1_filtered.sum() == 1
    assert df.is_cr2_filtered.sum() == 1
    assert df.is_call_rate_filtered.sum() == 2
    assert (~df.is_call_rate_filtered.astype("boolean")).sum() == 1
    # NOTE: Need to set type as "boolean" (instead of `bool`) to be able to
    # invert the boolean array containing missing values.
