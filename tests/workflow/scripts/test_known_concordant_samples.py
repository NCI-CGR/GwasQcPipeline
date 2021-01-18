import pandas as pd
import pytest
from typer.testing import CliRunner

from cgr_gwas_qc.models.config import Config
from cgr_gwas_qc.testing.data import RealData
from cgr_gwas_qc.workflow.scripts.known_concordant_samples import (
    app,
    read_pairwise_sample_concordance,
)

runner = CliRunner()


@pytest.mark.real_data
def test_unknown_concordant_samples(tmp_path, monkeypatch):
    """Make sure an unknown can be saved out.

    I am forcing 1 unknown pair of concordant samples. Here I am just testing
    that they are indeed written out.
    """
    # GIVEN: Real test data and the concordance threshold
    data_repo = RealData()
    threshold = str(Config().software_params.dup_concordance_cutoff)

    # Wrap the function that calculates concordance and replace the pair of
    # samples with the least concordance with a concordance value of 0.99. In
    # other words, force one unknown concordant sample.
    def mock_get_concordance(*args, **kwargs):
        df = read_pairwise_sample_concordance(*args, **kwargs)
        df.loc[df.concordance.argmin(), "concordance"] = 0.99
        return df

    # Monkeypatch the function in so when we run the script it uses our wrapped version.
    monkeypatch.setattr(
        "cgr_gwas_qc.workflow.scripts.known_concordant_samples.read_pairwise_sample_concordance",
        mock_get_concordance,
    )

    # WHEN: I run the script with legacy inputs
    results = runner.invoke(
        app,
        [
            (data_repo / "original_data/manifest_full.csv").as_posix(),
            (
                data_repo / "production_outputs/plink_filter_call_rate_2/samples_filter2.imiss"
            ).as_posix(),
            (data_repo / "production_outputs/ibd/samples.genome").as_posix(),
            (tmp_path / "known.csv").as_posix(),
            (tmp_path / "known_qc.csv").as_posix(),
            (tmp_path / "known_study.csv").as_posix(),
            (tmp_path / "unknown.csv").as_posix(),
            threshold,
        ],
        catch_exceptions=False,
    )
    assert results.exit_code == 0  # Make sure it ran successfully

    # THEN: There should be 1 unknown concordant sample.
    obs_unknown = pd.read_csv(tmp_path / "unknown.csv")
    assert obs_unknown.shape[0] == 1
