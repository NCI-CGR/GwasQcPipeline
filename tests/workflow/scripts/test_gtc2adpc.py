import numpy as np
import pytest
from numpy import isclose

from cgr_gwas_qc.parsers.illumina import AdpcReader
from cgr_gwas_qc.workflow.scripts import gtc2adpc


@pytest.mark.regression
def test_gtc2adpc_fake_data(fake_data_cache, tmp_path):
    # GIVEN: A fake GTC and BPM file
    obs_adpc = tmp_path / "test.adpc.bin"
    gtc2adpc.main(
        fake_data_cache / "illumina/gtc/small_genotype.gtc",
        fake_data_cache / "illumina/bpm/small_manifest.bpm",
        obs_adpc,
    )

    # THEN: The adpc.bin file is almost identical are about equal. The output
    # contains floats so we need to do fuzzy matching.
    expected_adpc = fake_data_cache / "small.adpc.bin"
    with AdpcReader(expected_adpc) as expected, AdpcReader(obs_adpc) as obs:
        for exp_row, obs_row in zip(expected, obs):
            # Integers should be the exactly the same
            assert exp_row.x_raw == obs_row.x_raw
            assert exp_row.y_raw == obs_row.y_raw
            assert exp_row.genotype == obs_row.genotype

            # Floats should be really close
            assert isclose(exp_row.x_norm, obs_row.x_norm)
            assert isclose(exp_row.y_norm, obs_row.y_norm)
            assert isclose(exp_row.genotype_score, obs_row.genotype_score)


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.slow
def test_gtc2adpc_issue_31(real_data_cache, tmp_path):
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
    obs_adpc = tmp_path / "test.adpc.bin"
    gtc2adpc.main(
        real_data_cache / "original_data/202917110153_R02C02.gtc",
        real_data_cache / "reference_data/GSAMD-24v1-0_20011747_A1.bpm",
        obs_adpc,
    )

    # THEN:
    expected_adpc = real_data_cache / "legacy_outputs/contam/SB034327_PC26469_C09.adpc.bin"
    with AdpcReader(expected_adpc) as expected, AdpcReader(obs_adpc) as obs:
        na_count = 0
        for exp_row, obs_row in zip(expected, obs):
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
                assert isclose(exp_row.x_norm, obs_row.x_norm)
                assert isclose(exp_row.y_norm, obs_row.y_norm)

                # genotypes scores are very close
                assert isclose(exp_row.genotype_score, obs_row.genotype_score)

        # The number of missing for this file is 236
        assert na_count == 236
