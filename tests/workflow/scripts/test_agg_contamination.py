from textwrap import dedent

import pytest

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
def agg_df(contam_csv, intensity_csv, imiss_file):
    return agg_contamination.build(contam_csv, intensity_csv, imiss_file)


def test_mask_low_intensity(agg_df, software_params):
    """Check that two NAs are added.

    SP00003 should be NA b/c of low intensity.
    SP00004 should be NA b/c of missing in imiss file.
    """
    obs_df = agg_contamination._mask_low_intensity(agg_df, software_params.intensity_threshold)
    assert obs_df["%Mix"].isna().sum() == 2


def test_flag_contaminated(agg_df, software_params):
    """Check that one sample is contaminated.

    SP00002 should be contaminated.
    """
    obs_df = agg_contamination._flag_contaminated(agg_df, software_params.contam_threshold)
    assert 1 == obs_df["is_contaminated"].sum()
    assert "SP00002" == obs_df[obs_df.is_contaminated].squeeze().Sample_ID
