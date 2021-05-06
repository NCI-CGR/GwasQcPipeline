import pytest
from pytest_mock import MockerFixture

from cgr_gwas_qc.testing import chdir
from cgr_gwas_qc.testing.data import FakeData


@pytest.mark.parametrize(
    "cgems,biowulf,cluster_profile",
    [(True, False, None), (False, True, None), (False, False, "./test")],
)
def test_check_exclusive_options_no_error(cgems, biowulf, cluster_profile):
    from cgr_gwas_qc.cli.submit import check_exclusive_options

    check_exclusive_options(cgems, biowulf, cluster_profile)


@pytest.mark.parametrize(
    "cgems,biowulf,cluster_profile",
    [(True, True, "./test"), (True, True, None), (True, False, "./test"), (False, True, "./test")],
)
def test_check_exclusive_options_raises_error(cgems, biowulf, cluster_profile):
    from click.exceptions import Exit as ClickExit

    from cgr_gwas_qc.cli.submit import check_exclusive_options

    with pytest.raises(ClickExit):
        check_exclusive_options(cgems, biowulf, cluster_profile)


@pytest.mark.parametrize("cluster", ["cgems", pytest.param("biowulf", marks=pytest.mark.xfail)])
def test_get_profile(cluster):
    from pathlib import Path

    from cgr_gwas_qc.cli.submit import get_profile

    profile = Path(get_profile(cluster))
    assert profile.exists() & profile.is_dir()


def test_check_custom_cluster_profile(tmp_path):
    from cgr_gwas_qc.cli.submit import check_custom_cluster_profile

    cluster_profile = tmp_path / "test_present"
    cluster_profile.mkdir()
    queue = "all"
    submission_cmd = "qsub"
    assert cluster_profile.resolve().as_posix() == check_custom_cluster_profile(
        cluster_profile, queue, submission_cmd
    )


def test_check_custom_cluster_profile_no_profile(tmp_path):
    from cgr_gwas_qc.cli.submit import check_custom_cluster_profile

    # Missing profile directory
    with pytest.raises(ValueError):
        cluster_profile = tmp_path / "test_present"
        queue, cmd = "all", "qsub"
        check_custom_cluster_profile(cluster_profile, queue, cmd)


def test_check_custom_cluster_profile_no_queue(tmp_path):
    from cgr_gwas_qc.cli.submit import check_custom_cluster_profile

    # Missing profile directory
    with pytest.raises(ValueError):
        cluster_profile = tmp_path / "test_present"
        cluster_profile.mkdir()
        queue, cmd = None, "qsub"
        check_custom_cluster_profile(cluster_profile, queue, cmd)


def test_check_custom_cluster_profile_no_cmd(tmp_path):
    from cgr_gwas_qc.cli.submit import check_custom_cluster_profile

    # Missing profile directory
    with pytest.raises(ValueError):
        cluster_profile = tmp_path / "test_present"
        cluster_profile.mkdir()
        queue, cmd = "all", None
        check_custom_cluster_profile(cluster_profile, queue, cmd)


@pytest.mark.parametrize("sample_size,group_size", [(10, 500), (101, 500), (10_000, 1000)])
def test_get_group_size_mocks(sample_size, group_size):
    from cgr_gwas_qc.cli.submit import get_group_size

    assert group_size == get_group_size(sample_size)


def test_get_per_sample_rules():
    from cgr_gwas_qc.cli.submit import get_per_sample_rules

    per_sample_rules = get_per_sample_rules()
    assert all(x.startswith("per_sample_") for x in per_sample_rules)
    assert 4 == len(per_sample_rules)


@pytest.mark.parametrize("sample_size,group_size", [(10, 500), (101, 500), (10_000, 1000)])
def test_get_grouping_settings(sample_size, group_size):
    from cgr_gwas_qc.cli.submit import get_grouping_settings

    expected = (
        "--group-components "
        f"per_sample_gtc_to_adpc={group_size} "
        f"per_sample_gtc_to_ped={group_size} "
        f"per_sample_median_idat_intensity={group_size} "
        f"per_sample_verifyIDintensity_contamination={group_size}"
    )
    assert expected == get_grouping_settings(sample_size)


def test_create_submission_script_cgems(tmp_path):
    import os

    from cgr_gwas_qc.cli.submit import create_submission_script

    with chdir(tmp_path):
        payload = {
            "python_executable": "python",
            "working_dir": os.getcwd(),
            "cgems": True,
            "biowulf": False,
            "time_h": 12,
            "queue": "all.q",
            "profile": "test_profile",
            "local_tasks": 1,
            "local_mem_mb": 500,
            "group_options": "",
        }
        create_submission_script(payload)

    assert (
        "#$ -N GwasQcPipeline" in (tmp_path / ".snakemake/GwasQcPipeline_submission.sh").read_text()
    )


def test_create_submission_script_biowulf(tmp_path):
    import os

    from cgr_gwas_qc.cli.submit import create_submission_script

    with chdir(tmp_path):
        payload = {
            "python_executable": "python",
            "working_dir": os.getcwd(),
            "cgems": False,
            "biowulf": True,
            "time_h": 12,
            "queue": "all.q",
            "profile": "test_profile",
            "group_options": "",
        }
        create_submission_script(payload)

    assert (
        '#SBATCH --job-name="GwasQcPipeline"'
        in (tmp_path / ".snakemake/GwasQcPipeline_submission.sh").read_text()
    )


@pytest.mark.parametrize(
    "cluster,cmd", [("cgems", "qsub"), ("biowulf", "sbatch"), ("custom", "pbs")]
)
def test_run_submit_with_right_command(cluster, cmd, tmp_path, mocker: MockerFixture):
    from cgr_gwas_qc.cli import submit

    if cluster == "cgems":
        cgems, biowulf, cluster_profile, queue, submission_cmd = True, False, None, None, None
    elif cluster == "biowulf":  # Biowulf
        cgems, biowulf, cluster_profile, queue, submission_cmd = False, True, None, None, None
    else:
        profile_dir = tmp_path / "test_profile"
        profile_dir.mkdir()
        cgems, biowulf, cluster_profile, queue, submission_cmd = (
            False,
            False,
            profile_dir,
            "all",
            cmd,
        )

    spy = mocker.patch("cgr_gwas_qc.cli.submit.sp.check_output")
    FakeData(tmp_path).make_cgr_sample_sheet().make_config()
    with chdir(tmp_path):
        submit.main(
            cgems=cgems,
            biowulf=biowulf,
            cluster_profile=cluster_profile,
            queue=queue,
            submission_cmd=submission_cmd,
            time_hr=12,
            dry_run=False,
        )

    spy.assert_called_once_with([cmd, ".snakemake/GwasQcPipeline_submission.sh"])
