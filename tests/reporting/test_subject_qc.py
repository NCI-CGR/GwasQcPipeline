import shutil
from io import StringIO
from pathlib import Path
from typing import Tuple

import pandas as pd
import pytest

# real_config, sample_sheet_full, sample_qc_df, population_qc_df


@pytest.mark.real_data
def test_UnExpectedReplicates(real_data_cache, subject_qc_df):
    from cgr_gwas_qc.reporting.subject_qc import UnExpectedReplicates

    unexpected_replicates_df = pd.read_csv(
        real_data_cache / "dev_outputs/sample_level/concordance/UnknownReplicates.csv"
    )
    unexpected = UnExpectedReplicates.construct(subject_qc_df, unexpected_replicates_df)

    assert 0 == unexpected.num_unexpected_replicates


@pytest.mark.real_data
def test_SexVerification(real_cfg, subject_qc_df):
    from cgr_gwas_qc.reporting.subject_qc import SexVerification

    sv = SexVerification.construct(real_cfg.ss, subject_qc_df, Path("test.png"))

    assert 3 == sv.num_sex_discordant

    # The markdown table should have the same number of rows
    obs_df = (
        pd.read_csv(StringIO(sv.table), sep="|", skipinitialspace=True)
        .dropna(axis=1, how="all")
        .iloc[1:, :]
    )
    obs_df.columns = obs_df.columns.str.strip()
    assert 3 == obs_df.shape[0]


@pytest.mark.real_data
def test_Relatedness(agg_population_qc_df):
    from cgr_gwas_qc.reporting.subject_qc import Relatedness

    rl = Relatedness.construct(agg_population_qc_df)

    assert 0 == rl.num_related_subjects
    assert 0 == rl.num_qc_families


@pytest.fixture(params=[1, 2, 5])
def num_and_png_dir(request, tmp_path, fake_image) -> Tuple[int, Path]:
    for i in range(request.param):
        shutil.copy(fake_image, tmp_path / f"Pop{i}.png")

    return request.param, tmp_path


@pytest.mark.real_data
def test_Autosomal(agg_population_qc_df, num_and_png_dir):
    from cgr_gwas_qc.reporting.subject_qc import Autosomal

    num_pops, png_dir = num_and_png_dir
    az = Autosomal.construct(agg_population_qc_df, png_dir)

    assert 1 == az.num_populations_analyzed
    assert 0 == az.num_subjects_excluded
    assert num_pops == len(az.panels)


def test_Pca(num_and_png_dir):
    from cgr_gwas_qc.reporting.subject_qc import Pca

    num_pops, png_dir = num_and_png_dir
    pca = Pca.construct(png_dir)

    assert num_pops == len(pca.panels)


def test_Hwe(num_and_png_dir):
    from cgr_gwas_qc.reporting.subject_qc import Hwe

    num_pops, png_dir = num_and_png_dir
    hwe = Hwe.construct(png_dir)

    assert num_pops == len(hwe.panels)
