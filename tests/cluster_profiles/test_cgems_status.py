import pytest

from cgr_gwas_qc.cluster_profiles.cgems import cgems_status


@pytest.mark.parametrize("job_id,expected_status", [(2, "running"), (4, "failed"), (22, None)])
def test_check_queue(qsub, job_id, expected_status):
    assert expected_status == cgems_status.check_queue(job_id)


@pytest.mark.parametrize("job_id,expected_status", [(1, "success"), (0, "failed")])
def test_check_job_history(qsub, job_id, expected_status):
    assert expected_status == cgems_status.check_job_history(job_id)
