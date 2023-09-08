import os
from pathlib import Path
from shutil import which

import pytest
import snakemake


@pytest.mark.skipif(which("qsub") is None, reason="Not an SGE system with qsub command to test")
@pytest.mark.skipif(
    os.environ.get("GITHUB_ACTIONS", "false") == "true", reason="Broken on GitHub Actions"
)
def test_cgems_submission_end_to_end(tmp_path, qsub):
    # tmp_path = "src/cgr_gwas_qc/testing"
    print("!!!!!!!!!")
    print(tmp_path)
    print("!!!!!!!!!")
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
