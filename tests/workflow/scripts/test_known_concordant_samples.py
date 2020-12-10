import pandas as pd
from numpy import isclose
from typer.testing import CliRunner

from cgr_gwas_qc.models.config import Config
from cgr_gwas_qc.testing.data import RealData
from cgr_gwas_qc.workflow.scripts.known_concordant_samples import (
    app,
    read_pairwise_sample_concordance,
)

runner = CliRunner()


def test_known_concordant_samples(tmp_path):
    # GIVEN: Real data
    data_repo = RealData()
    cfg = Config()

    # WHEN: Run the known_concordan_samples.py with real data inputs
    results = runner.invoke(
        app,
        [
            (data_repo / "original_data/manifest_full.csv").as_posix(),
            (
                data_repo / "production_outputs/plink_filter_call_rate_2/samples_filter2.imiss"
            ).as_posix(),
            (data_repo / "production_outputs/ibd/samples.genome").as_posix(),
            str(cfg.software_params.dup_concordance_cutoff),
            (tmp_path / "known.csv").as_posix(),
            (tmp_path / "unknown.csv").as_posix(),
        ],
        catch_exceptions=False,
    )

    # THEN: The script should run without error
    assert results.exit_code == 0

    # The values of the known replicates should have the same (concordance and PI_HAT)
    obs_known = pd.read_csv(tmp_path / "known.csv")
    exp_known = pd.read_csv(data_repo / "production_outputs/concordance/KnownReplicates.csv")
    for obs_r in obs_known.itertuples():
        exp_r = exp_known.query("Sample_ID1 == @obs_r.Sample_ID1 & Sample_ID2 == @obs_r.Sample_ID2")
        assert isclose(obs_r.concordance, exp_r.Concordance.values[0])
        assert isclose(obs_r.PI_HAT, exp_r.PI_HAT.values[0])

    # There were no unknown replicates in the real data.
    obs_unknown = pd.read_csv(tmp_path / "unknown.csv")
    assert obs_unknown.shape[0] == 0


def test_unknown_concordant_samples(tmp_path, monkeypatch):
    # GIVEN: Real data and 1 unknown concordant pair of samples.
    data_repo = RealData()
    cfg = Config()

    def mock_get_concordance(*args, **kwargs):
        df = read_pairwise_sample_concordance(*args, **kwargs)
        df.loc[df.concordance.argmin(), "concordance"] = 0.99
        return df

    monkeypatch.setattr(
        "cgr_gwas_qc.workflow.scripts.known_concordant_samples.read_pairwise_sample_concordance",
        mock_get_concordance,
    )

    # WHEN: Run the known_concordan_samples.py with real data inputs
    results = runner.invoke(
        app,
        [
            (data_repo / "original_data/manifest_full.csv").as_posix(),
            (
                data_repo / "production_outputs/plink_filter_call_rate_2/samples_filter2.imiss"
            ).as_posix(),
            (data_repo / "production_outputs/ibd/samples.genome").as_posix(),
            str(cfg.software_params.dup_concordance_cutoff),
            (tmp_path / "known.csv").as_posix(),
            (tmp_path / "unknown.csv").as_posix(),
        ],
        catch_exceptions=False,
    )

    # THEN: The script should run without error
    assert results.exit_code == 0

    # There should be 1 unknown concordant sample.
    obs_unknown = pd.read_csv(tmp_path / "unknown.csv")
    assert obs_unknown.shape[0] == 1
