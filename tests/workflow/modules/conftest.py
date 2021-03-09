import pytest

from cgr_gwas_qc.testing.data import FakeData, RealData


@pytest.fixture(params=["fake", pytest.param("real", marks=pytest.mark.real_data)])
def gtc_data_repo(request, tmp_path, conda_envs):
    """Creates a DataRepo for Fake and Real data.

    The DataRepo working_dir contains:
        - all conda environments
        - sample sheet
        - user files (GTC, IDATs; fake data only)
        - config.yml contains:
            - reference files (BPM, VCF, TBI)
            - user files (gtc_pattern, idat_pattern)
            - software params (defaults)
            - workflow params (defaults)

    Returns:
        A FakeData or RealData data cache.
    """
    conda_envs.copy_all_envs(tmp_path)

    if request.param == "real":
        return (
            RealData(tmp_path)
            .add_sample_sheet(full_sample_sheet=False)
            .add_user_files(copy=False, entry_point="gtc")
            .make_config()
        )
    else:
        return FakeData(tmp_path).add_sample_sheet().add_user_files(entry_point="gtc").make_config()
