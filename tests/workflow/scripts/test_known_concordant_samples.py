import pandas as pd
import pytest
from typer.testing import CliRunner

from cgr_gwas_qc import load_config
from cgr_gwas_qc.testing import chdir
from cgr_gwas_qc.testing.data import RealData
from cgr_gwas_qc.workflow.scripts.known_concordant_samples import app

runner = CliRunner()


@pytest.mark.real_data
def test_unknown_concordant_samples(tmp_path):
    """Make sure an unknown can be saved out.

    I am forcing 1 unknown pair of concordant samples. Here I am just testing
    that they are indeed written out.
    """
    # GIVEN: Real test data and the concordance threshold
    data_cache = RealData(tmp_path).add_sample_sheet().make_config()

    with chdir(tmp_path):
        cfg = load_config()

    threshold = str(cfg.config.software_params.dup_concordance_cutoff)

    # Create a concordance table with 1 unknown concordant sample
    (
        pd.read_csv(data_cache / "production_outputs/ibd/samples.genome", delim_whitespace=True)
        .assign(concordance=lambda x: x.IBS2 / (x.IBS0 + x.IBS1 + x.IBS2))
        .assign(  # Replace the minimum concordance with .99 to cause an unexpected replicate
            concordance=lambda x: x.concordance.map(
                lambda y: 0.99 if y == x.concordance.min() else y
            )
        )
        .reindex(["IID1", "IID2", "PI_HAT", "concordance"], axis=1)
        .to_csv(tmp_path / "one_unknown_concordant.csv")
    )

    # WHEN: I run the script with legacy inputs
    results = runner.invoke(
        app,
        [
            (data_cache / "original_data/manifest_full.csv").as_posix(),
            (
                data_cache / "production_outputs/plink_filter_call_rate_2/samples_filter2.imiss"
            ).as_posix(),
            (tmp_path / "one_unknown_concordant.csv").as_posix(),
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
