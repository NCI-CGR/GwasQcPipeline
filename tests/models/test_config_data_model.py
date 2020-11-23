from pathlib import Path

import pytest
from pydantic import ValidationError

from cgr_gwas_qc.models.config import (
    Config,
    ReferenceFiles,
    SoftwareParams,
    UserFiles,
    WorkflowParams,
)
from cgr_gwas_qc.testing import chdir
from cgr_gwas_qc.version import __version__


def test_reference_files_must_exist():
    with pytest.raises(ValidationError):
        ReferenceFiles(
            illumina_manifest_file="test.bpm",
            thousand_genome_vcf="test.vcf.gz",
            thousand_genome_tbi="test.vcf.gz.tbi",
        )


def test_reference_files(tmp_path):
    with chdir(tmp_path):
        [Path(x).touch() for x in ("test.bpm", "test.vcf.gz", "test.vcf.gz.tbi")]

        refs = ReferenceFiles(
            illumina_manifest_file="test.bpm",
            thousand_genome_vcf="test.vcf.gz",
            thousand_genome_tbi="test.vcf.gz.tbi",
        )

        assert refs.illumina_manifest_file == Path("test.bpm")
        assert refs.thousand_genome_vcf == Path("test.vcf.gz")
        assert refs.thousand_genome_tbi == Path("test.vcf.gz.tbi")


def test_user_files_defaults():
    user_files = UserFiles()
    assert user_files.gtc_pattern is None
    assert user_files.idat_pattern is None
    assert user_files.ped is None
    assert user_files.map is None
    assert user_files.bed is None
    assert user_files.bim is None
    assert user_files.fam is None


def test_software_params_defaults():
    software_params = SoftwareParams()
    assert software_params.samp_cr_1 == 0.8
    assert software_params.snp_cr_1 == 0.8
    assert software_params.samp_cr_2 == 0.95
    assert software_params.snp_cr_2 == 0.95
    assert software_params.ld_prune_r2 == 0.1
    assert software_params.maf_for_ibd == 0.2
    assert software_params.ibd_pi_hat_cutoff == 0.95
    assert software_params.dup_concordance_cutoff == 0.95
    assert software_params.contam_threshold == 0.2
    assert software_params.contam_population == "AF"
    assert software_params.pi_hat_threshold == 0.2
    assert software_params.autosomal_het_threshold == 0.1
    assert software_params.minimum_pop_subjects == 50
    assert software_params.control_hwp_threshold == 50
    assert software_params.strand == "top"


def test_workflow_params_defaults():
    workflow_params = WorkflowParams()
    assert workflow_params.expected_sex_col_name == "Expected_Sex"
    assert workflow_params.remove_contam is True
    assert workflow_params.remove_sex_discordant is True
    assert workflow_params.remove_rep_discordant is True
    assert workflow_params.remove_unexpected_rep is True


def test_basic_config(tmp_path):
    with chdir(tmp_path):
        [
            Path(x).touch()
            for x in ("test.bpm", "test.vcf.gz", "test.vcf.gz.tbi", "sample_sheet.csv")
        ]
        cfg = Config(
            project_name="Test Project",
            sample_sheet="sample_sheet.csv",
            reference_files=ReferenceFiles(
                illumina_manifest_file="test.bpm",
                thousand_genome_vcf="test.vcf.gz",
                thousand_genome_tbi="test.vcf.gz.tbi",
            ),
            user_files=UserFiles(),
            software_params=SoftwareParams(),
            workflow_params=WorkflowParams(),
        )

        assert cfg.project_name == "Test Project"
        assert cfg.sample_sheet == Path("sample_sheet.csv")
        assert cfg.pipeline_version == __version__


@pytest.mark.parametrize("value", [-0.01, 0, 1.01])
def test_call_rate_param_limits(value):
    """Make sure call rate validation is working.

    Again more a sanity check that pydantic is behaving as expected. Sample
    and SNP call rate filters should between 0-1. Values ≤0 or >1 should
    raise a ``pydantic.ValidationError``.
    """

    with pytest.raises(ValidationError):
        SoftwareParams(samp_cr_1=value)

    with pytest.raises(ValidationError):
        SoftwareParams(snp_cr_1=value)

    with pytest.raises(ValidationError):
        SoftwareParams(samp_cr_2=value)

    with pytest.raises(ValidationError):
        SoftwareParams(snp_cr_2=value)
