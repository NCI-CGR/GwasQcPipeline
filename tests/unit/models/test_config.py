from pathlib import Path

import pytest
from pydantic import ValidationError

from cgr_gwas_qc import load_config
from cgr_gwas_qc.models import config as cfg_models
from cgr_gwas_qc.testing import chdir
from cgr_gwas_qc.version import __version__


def test_load_config_from_yml(gtc_working_dir):
    """Make sure loading YAML config gives expected results.

    This is just a sanity check. Assuming pydantic works as expected, then
    this should always pass unless we change default values.
    """
    with chdir(gtc_working_dir):
        cfg = load_config()
        assert cfg.config.pipeline_version == __version__
        assert cfg.config.project_name == "Test Project"
        assert cfg.config.sample_sheet == Path("example_sample_sheet.csv")

        # Reference Files
        assert cfg.config.reference_files.illumina_manifest_file == Path("small_manifest.bpm")
        assert cfg.config.reference_files.thousand_genome_vcf == Path("small_1KG.vcf.gz")
        assert cfg.config.reference_files.thousand_genome_tbi == Path("small_1KG.vcf.gz.tbi")

        # User Files
        # Since we are using the ``gtc_working_dir`` only GTC and IDAT should be defined.
        assert cfg.config.user_files.gtc_pattern == "{Sample_ID}.gtc"
        assert cfg.config.user_files.idat_pattern.red == "{Sample_ID}_Red.idat"
        assert cfg.config.user_files.idat_pattern.green == "{Sample_ID}_Grn.idat"

        assert cfg.config.user_files.ped is None
        assert cfg.config.user_files.map is None

        assert cfg.config.user_files.bed is None
        assert cfg.config.user_files.bim is None
        assert cfg.config.user_files.fam is None

        # Software Params
        assert cfg.config.software_params.samp_cr_1 == 0.8
        assert cfg.config.software_params.snp_cr_1 == 0.8
        assert cfg.config.software_params.samp_cr_2 == 0.95
        assert cfg.config.software_params.snp_cr_2 == 0.95
        assert cfg.config.software_params.ld_prune_r2 == 0.1
        assert cfg.config.software_params.maf_for_ibd == 0.2
        assert cfg.config.software_params.ibd_pi_hat_cutoff == 0.95
        assert cfg.config.software_params.dup_concordance_cutoff == 0.95
        assert cfg.config.software_params.contam_threshold == 0.2
        assert cfg.config.software_params.contam_population == "AF"
        assert cfg.config.software_params.pi_hat_threshold == 0.2
        assert cfg.config.software_params.autosomal_het_threshold == 0.1
        assert cfg.config.software_params.minimum_pop_subjects == 50
        assert cfg.config.software_params.control_hwp_threshold == 50
        assert (
            cfg.config.software_params.strand == "top"
        )  # Note we convert strand to lowercase on import

        # Workflow Params
        assert cfg.config.workflow_params.expected_sex_col_name == "Expected_Sex"
        assert cfg.config.workflow_params.remove_contam is True
        assert cfg.config.workflow_params.remove_sex_discordant is True
        assert cfg.config.workflow_params.remove_rep_discordant is True
        assert cfg.config.workflow_params.remove_unexpected_rep is True


@pytest.mark.parametrize("value", [-0.01, 0, 1.01])
def test_call_rate_param_limits(value):
    """Make sure call rate validation is working.

    Again more a sanity check that pydantic is behaving as expected. Sample
    and SNP call rate filters should between 0-1. Values ≤0 or >1 should
    raise a ``pydantic.ValidationError``.
    """

    with pytest.raises(ValidationError):
        cfg_models.SoftwareParams(samp_cr_1=value)

    with pytest.raises(ValidationError):
        cfg_models.SoftwareParams(snp_cr_1=value)

    with pytest.raises(ValidationError):
        cfg_models.SoftwareParams(samp_cr_2=value)

    with pytest.raises(ValidationError):
        cfg_models.SoftwareParams(snp_cr_2=value)
