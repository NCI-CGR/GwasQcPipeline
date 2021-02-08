import pytest
from pytest_mock import MockerFixture

from cgr_gwas_qc.testing import chdir


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


@pytest.mark.parametrize("cluster", ["cgems", "biowulf"])
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

    spy = mocker.patch("cgr_gwas_qc.cli.submit.sp.run")
    with chdir(tmp_path):
        submit.main(
            cgems=cgems,
            biowulf=biowulf,
            cluster_profile=cluster_profile,
            queue=queue,
            submission_cmd=submission_cmd,
            time_h=12,
        )

    spy.assert_called_once_with(
        [cmd, ".snakemake/GwasQcPipeline_submission.sh"], shell=True, check=True
    )
