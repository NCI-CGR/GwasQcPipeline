from pathlib import Path

from cgr_gwas_qc.models.config import Config
from cgr_gwas_qc.testing import chdir
from cgr_gwas_qc.yaml import load


def test_config_model(gtc_working_dir):
    with chdir(gtc_working_dir):
        cfg = Config(**load(Path("config.yml")))
        assert cfg.project_name == "Test Project"
        assert cfg.sample_sheet == Path("example_sample_sheet.csv")
