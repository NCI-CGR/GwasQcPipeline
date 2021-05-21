from textwrap import dedent

import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal, assert_series_equal

from cgr_gwas_qc.workflow.scripts import agg_contamination


@pytest.fixture
def contam_csv(tmp_path):
    outfile = tmp_path / "contam.csv"
    outfile.write_text(
        dedent(
            """\
            Sample_ID,%Mix,LLK,LKK0
            SP00001,0.01,-2,0.2
            SP00002,0.20,-2,0.2
            SP00003,0.01,-2,0.2
            SP00004,0.01,-2,0.2
            SP00005,0.01,-2,0.2
            """
        )
    )

    return outfile


@pytest.fixture
def intensity_csv(tmp_path):
    outfile = tmp_path / "intensity.csv"
    outfile.write_text(
        dedent(
            """\
            Sample_ID,Chip_ID,median_intensity
            SP00001,test1,7000.0
            SP00002,test1,7000.0
            SP00003,test1,4000.0
            SP00004,test1,7000.0
            SP00005,test1,4000.0
            """
        )
    )
    return outfile


@pytest.fixture
def imiss_file(tmp_path):
    outfile = tmp_path / "test.imiss"
    outfile.write_text(
        dedent(
            """\
            FID   IID   MISS_PHENO   N_MISS   N_GENO   F_MISS
            SB00001  SP00001  Y  7000  700000  0.01
            SB00002  SP00002  Y  7000  700000  0.01
            SB00003  SP00003  Y  7000  700000  0.01
            """
        )
    )
    return outfile


@pytest.fixture
def contam_df():
    return pd.DataFrame(
        [
            #   0        1     2   3      4        5     6    7      8       9     10
            ("SP00001", 0.01, -2, 0.2, "test1", 7000.0, "Y", 7000, 700000, 0.01, False),
            ("SP00002", 0.20, -2, 0.2, "test1", 7000.0, "Y", 7000, 700000, 0.01, True),
            ("SP00003", 0.01, -2, 0.2, "test1", 4000.0, "Y", 7000, 700000, 0.01, False),
            ("SP00004", 0.01, -2, 0.2, "test1", 7000.0, np.nan, np.nan, np.nan, np.nan, False),
            ("SP00005", 0.01, -2, 0.2, "test1", 4000.0, np.nan, np.nan, np.nan, np.nan, False),
        ],
        columns=(
            "Sample_ID",  # 0
            "%Mix",  # 1
            "LLK",  # 2
            "LKK0",  # 3
            "Chip_ID",  # 4
            "median_intensity",  # 5
            "MISS_PHENO",  # 6
            "N_MISS",  # 7
            "N_GENO",  # 8
            "F_MISS",  # 9
            "is_contaminated",  # 10
        ),
    )


def test_agg_df(contam_csv, intensity_csv, imiss_file, contam_df):
    expected = contam_df.drop("is_contaminated", axis=1)
    obs = agg_contamination.build(contam_csv, intensity_csv, imiss_file)
    assert_frame_equal(expected, obs)


def test_mask_low_intensity(contam_df, software_params):
    """Check that two NAs are added.

    SP00003 should be NA b/c of low intensity.
    SP00004 should be NA b/c of missing in imiss file.
    """
    # GIVEN: contamination table
    obs_df = contam_df.drop("is_contaminated", axis=1)
    expected_df = contam_df.drop("is_contaminated", axis=1)

    # WHEN: I mask low intensity values
    threshold = software_params.intensity_threshold
    agg_contamination._mask_low_intensity(obs_df, threshold)

    # THEN: %MIX should be NA at low intensity
    expected_df.loc[4, "%Mix"] = pd.NA  # SP00005 b/c < threshold and missing F_MISS
    assert_frame_equal(expected_df, obs_df)


def test_flag_contaminated(contam_df, software_params):
    """Check that one sample is contaminated.

    SP00002 should be contaminated.
    """
    obs_df = (
        contam_df.drop("is_contaminated", axis=1)
        .pipe(agg_contamination._mask_low_intensity, software_params.intensity_threshold)
        .pipe(agg_contamination._flag_contaminated, software_params.contam_threshold)
    )

    assert_series_equal(contam_df.is_contaminated, obs_df.is_contaminated)


@pytest.mark.real_data
@pytest.mark.regression
def test_legacy_masking(real_data_cache, software_params, tmp_path):
    # GIVEN: legacy and dev workflow outputs
    legacy = pd.read_csv(
        real_data_cache / "legacy_outputs/all_contam/contam.csv",
        index_col=0,
        usecols=["ID", "%Mix"],
    ).squeeze()

    # WHEN: Run the agg script
    outfile = tmp_path / "test.csv"
    agg_contamination.main(
        real_data_cache / "dev_outputs/sample_level/contamination/verifyIDintensity.csv",
        real_data_cache / "dev_outputs/sample_level/contamination/median_idat_intensity.csv",
        real_data_cache / "dev_outputs/sample_level/call_rate_2/samples.imiss",
        software_params.intensity_threshold,
        software_params.contam_threshold,
        outfile,
    )
    dev = pd.read_csv(outfile, index_col=0, usecols=["Sample_ID", "%Mix"]).squeeze()

    # THEN: legacy and dev masking should be identical
    assert_series_equal(legacy.isna(), dev.isna(), check_names=False)
