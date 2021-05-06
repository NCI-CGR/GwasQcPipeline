import os
from pathlib import Path

import pytest
import snakemake


@pytest.mark.skipif(
    os.environ.get("GITHUB_ACTIONS", "false") == "true", reason="Broken on GitHub Actions"
)
def test_cgems_submission_end_to_end(tmp_path, qsub):
    with pytest.raises(SystemExit) as exc:
        args = [
            "-j",
            "1",
            "-s",
            "tests/data/job_scripts/basic.smk",
            "-d",
            tmp_path.as_posix(),
            "--profile",
            Path("src/cgr_gwas_qc/cluster_profiles/cgems").absolute().as_posix(),
        ]

        snakemake.main(args)

    assert exc.value.code == 0
