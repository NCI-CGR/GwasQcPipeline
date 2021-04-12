import pandas as pd
import pytest

from cgr_gwas_qc.config import ConfigMgr
from cgr_gwas_qc.testing.data import RealData
from cgr_gwas_qc.workflow.scripts.known_concordant_samples import main


@pytest.mark.real_data
def test_unknown_concordant_samples(real_cfg: ConfigMgr, tmp_path):
    """Make sure an unknown can be saved out.

    I am forcing 1 unknown pair of concordant samples. Here I am just testing
    that they are indeed written out.
    """
    # GIVEN: Real test data and the concordance threshold
    data_cache = RealData()

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
    main(
        real_cfg.root / "cgr_sample_sheet.csv",
        data_cache / "production_outputs/plink_filter_call_rate_2/samples_filter2.imiss",
        tmp_path / "one_unknown_concordant.csv",
        tmp_path / "known.csv",
        tmp_path / "known_qc.csv",
        tmp_path / "known_study.csv",
        tmp_path / "unknown.csv",
        real_cfg.config.software_params.dup_concordance_cutoff,
    )

    # THEN: There should be 1 unknown concordant sample.
    obs_unknown = pd.read_csv(tmp_path / "unknown.csv")
    assert obs_unknown.shape[0] == 1
