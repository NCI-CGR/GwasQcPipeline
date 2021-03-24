from textwrap import dedent

import numpy as np
import pytest

from cgr_gwas_qc.workflow.scripts.agg_population_qc_tables import add_metadata, aggregate_qc_tables


@pytest.fixture(scope="module")
def per_population_qc_done(tmp_path_factory):
    """Create a set of fake population level QC tables."""
    tmp_path = tmp_path_factory.mktemp("population_qc_tables")
    populations = {
        "European": ["U-141833-5", "R-023930-2", "Q-111481-8"],
        "African_American": ["D-099281-9", "X-126266-7", "R-043181-7"],
        "Asian": ["F-141988-8", "Y-041829-5", "Q-019204-7"],
    }
    population_filenames = []
    for pop_, subjects in populations.items():
        filename = tmp_path / f"{pop_}.csv"
        population_filenames.append(filename.as_posix())
        filename.write_text(
            dedent(
                f"""
                population,Subject_ID,data
                {pop_},{subjects[0]},{np.random.randint(0, 50)}
                {pop_},{subjects[1]},{np.random.randint(0, 50)}
                {pop_},{subjects[2]},{np.random.randint(0, 50)}
                """
            )
        )

    (tmp_path / "per_population_qc.done").write_text("\n".join(population_filenames))

    return tmp_path / "per_population_qc.done"


def test_aggregate_qc_tables(per_population_qc_done):
    """Test that fake population data is aggregated."""
    df = aggregate_qc_tables(per_population_qc_done)
    assert all(df.population.value_counts() == 3)
    assert df.Subject_ID.unique().shape[0] == 9


@pytest.mark.real_data
def test_add_metadata(sample_qc, per_population_qc_done):
    """Test adding sample level metadata"""
    df = aggregate_qc_tables(per_population_qc_done).pipe(add_metadata, filename=sample_qc)

    assert all(df.population.value_counts() == 3)
    assert df.Subject_ID.unique().shape[0] == 9
    assert 1 == df.case_control.isnull().sum()  # One subject was dropped from study so should be NA
    assert df.query("case_control == 'control'").shape[0] == 6
