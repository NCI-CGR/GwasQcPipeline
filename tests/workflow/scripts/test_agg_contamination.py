from textwrap import dedent

import pytest

from cgr_gwas_qc.workflow.scripts import agg_contamination


@pytest.fixture
def contam_files(tmp_path):
    data = dedent(
        """\
        #
        #
        ID  %Mix  LLK  LKK0
        -------------------
        0  {}  -2  0.2
        """
    )
    s1 = tmp_path / "SP00001.contam"
    s1.write_text(data.format(0.01))

    s2 = tmp_path / "SP00002.contam"
    s2.write_text(data.format(0.20))

    s3 = tmp_path / "SP00003.contam"
    s3.write_text(data.format(0.01))

    s4 = tmp_path / "SP00004.contam"
    s4.write_text(data.format(0.01))

    return s1, s2, s3


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
def agg_df(contam_files, intensity_csv, imiss_file):
    return agg_contamination.build(contam_files, intensity_csv, imiss_file)


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
    assert 1 == obs_df["is_ge_contam_threshold"].sum()
    assert "SP00002" == obs_df[obs_df.is_ge_contam_threshold].squeeze().Sample_ID
