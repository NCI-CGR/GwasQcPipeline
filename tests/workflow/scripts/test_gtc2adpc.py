import numpy as np
import pytest
from numpy import isclose
from typer.testing import CliRunner

from cgr_gwas_qc.parsers.illumina import AdpcReader
from cgr_gwas_qc.testing.data import FakeData, RealData
from cgr_gwas_qc.workflow.scripts.gtc2adpc import app

runner = CliRunner()


@pytest.fixture(scope="module")
def gtc_file():
    return FakeData._data_path / FakeData._test_gtc


@pytest.fixture(scope="module")
def bpm_file():
    return FakeData._data_path / FakeData._illumina_manifest_file


@pytest.mark.regression
def test_gtc2adpc(gtc_file, bpm_file, tmp_path):
    # GIVEN: A fake GTC and BPM file
    # WHEN: I run the script and get adpc.bin and snp_counts as outputs.
    results = runner.invoke(
        app,
        [gtc_file.as_posix(), bpm_file.as_posix(), (tmp_path / "adpc.bin").as_posix()],
    )
    assert results.exit_code == 0

    # THEN: The adpc.bin file is almost identical are about equal. The output
    # contains floats so we need to do fuzzy matching.
    obs_adpc = tmp_path / "adpc.bin"
    expected_adpc = "tests/data/small.adpc.bin"
    with AdpcReader(obs_adpc) as obs, AdpcReader(expected_adpc) as expect:
        for obs_row, exp_row in zip(obs, expect):
            # Integers should be the exactly the same
            assert obs_row.x_raw == exp_row.x_raw
            assert obs_row.y_raw == exp_row.y_raw
            assert obs_row.genotype == exp_row.genotype

            # Floats should be really close
            assert isclose(obs_row.x_norm, exp_row.x_norm)
            assert isclose(obs_row.y_norm, exp_row.y_norm)
            assert isclose(obs_row.genotype_score, exp_row.genotype_score)


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.slow
def test_gtc2adpc_issue_31(tmp_path):
    """Issue 31: Regression difference NA vs 0.0

    I found when using real data my version of the script returns `nan` and
    the production version of the script is returning ~0. After some digging
    I found some small differences in the Illumina parsing library when x or
    y raw are 0. The `nan` is the current Illumina implementation and is
    probably the right value to use.

    This test shows that when these do occur, they only occur when the raw
    values are 0, genotype score is 0, and the genotype call is 3 (Unknown).
    """
    # GIVEN: a real GTC and BPM file.
    data_repo = RealData()
    gtc_file = (data_repo / "original_data/202917110153_R02C02.gtc").as_posix()
    bpm_file = (data_repo / "reference_data/GSAMD-24v1-0_20011747_A1.bpm").as_posix()

    # WHEN: run the script
    results = runner.invoke(app, [gtc_file, bpm_file, (tmp_path / "test.adpc.bin").as_posix()])
    assert results.exit_code == 0

    # THEN:
    obs_adpc = tmp_path / "test.adpc.bin"
    expected_adpc = data_repo / "production_outputs/contam/SB034327_PC26469_C09.adpc.bin"

    with AdpcReader(obs_adpc) as obs, AdpcReader(expected_adpc) as expect:
        na_count = 0
        for obs_row, exp_row in zip(obs, expect):
            if np.isnan(obs_row.x_norm) & np.isnan(obs_row.y_norm):
                na_count += 1
                # if nan then the raw values should be 0 in obs and expected
                assert obs_row.x_raw == 0
                assert exp_row.x_raw == 0

                assert obs_row.y_raw == 0
                assert exp_row.y_raw == 0

                # if nan then the genotype score should be 0 in obs and expected
                assert obs_row.genotype_score == 0.0
                assert exp_row.genotype_score == 0.0

                # if nan then the genotype should be 3 (unknown) in obs and expected
                assert obs_row.genotype == 3
                assert exp_row.genotype == 3
            else:
                # otherwise the norm values are similar between the two scripts
                assert isclose(obs_row.x_norm, exp_row.x_norm)
                assert isclose(obs_row.y_norm, exp_row.y_norm)

                # genotypes scores are very close
                assert isclose(obs_row.genotype_score, exp_row.genotype_score)

        # The number of missing for this file is 236
        assert na_count == 236
