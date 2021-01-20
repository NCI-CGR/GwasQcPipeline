from pathlib import Path

import pandas as pd
import pytest

from cgr_gwas_qc.parsers.sample_sheet import SampleSheet
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
            .add_reference_files(copy=False)
            .add_user_files(copy=False, entry_point="gtc")
            .make_config()
        )
    else:
        return (
            FakeData(tmp_path)
            .add_sample_sheet()
            .add_reference_files(copy=False)
            .add_user_files(entry_point="gtc")
            .make_config()
        )


@pytest.mark.real_data
@pytest.fixture(scope="session")
def qc_summary(tmp_path_factory) -> Path:
    """Add new QC columns to legacy QC summary table.

    The legacy QC summary table is missing several columns that are not
    needed by downstream rules. This fixture adds the necessary columns.

    Returns:
        Path to an updated sample qc summary table (CSV).
    """
    from cgr_gwas_qc.workflow.scripts.sample_qc_report import (
        IDENTIFILER_FLAGS,
        QC_HEADER,
        _case_control_encoder,
        _find_study_subject_representative,
        _find_study_subject_with_no_representative,
        _identifiler_reason,
    )

    tmp_path = tmp_path_factory.mktemp("qc_table")
    data_cache = RealData()

    # Add Group By and Internal Control columns
    ss = (
        SampleSheet(data_cache / "original_data/manifest_full.csv")
        .add_group_by_column("PI_Subject_ID")
        .data.assign(Internal_Control=lambda x: x.Sample_Group == "sVALD-001")
        .reindex(
            ["Sample_ID", "Group_By_Subject_ID", "Internal_Control", "Case/Control_Status"], axis=1
        )
    )

    legacy_qc_table = pd.read_csv(data_cache / "production_outputs/all_sample_qc.csv").merge(
        ss, on="Sample_ID", how="left"
    )

    # Use functions from QC script to add other columns
    legacy_qc_table["Case/Control_Status"] = legacy_qc_table["Case/Control_Status"].map(
        _case_control_encoder
    )
    legacy_qc_table["Identifiler_Reason"] = _identifiler_reason(legacy_qc_table, IDENTIFILER_FLAGS)
    legacy_qc_table["Subject_Representative"] = _find_study_subject_representative(legacy_qc_table)
    legacy_qc_table["Subject_Dropped_From_Study"] = _find_study_subject_with_no_representative(
        legacy_qc_table
    )

    legacy_qc_table.reindex(QC_HEADER, axis=1).to_csv(tmp_path / "qc.csv", index=False)
    return tmp_path / "qc.csv"
