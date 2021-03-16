from pathlib import Path

import pytest
from pydantic import ValidationError

from cgr_gwas_qc.models.config import (
    Idat,
    ReferenceFiles,
    SoftwareParams,
    UserFiles,
    WorkflowParams,
)


def test_reference_files():
    # GIVEN-WHEN: we load the reference files model with some values
    refs = ReferenceFiles(
        illumina_manifest_file="test.bpm",
        thousand_genome_vcf="test.vcf.gz",
        thousand_genome_tbi="test.vcf.gz.tbi",
    )

    # THEN: those values are accessible and convert to paths.
    assert refs.illumina_manifest_file == Path("test.bpm")
    assert refs.thousand_genome_vcf == Path("test.vcf.gz")
    assert refs.thousand_genome_tbi == Path("test.vcf.gz.tbi")


def test_user_files_defaults():
    # GIVEN-WHEN: we load the user files model without any inputs
    user_files = UserFiles()

    # THEN: all values default to None
    assert user_files.output_pattern == "{prefix}/{file_type}.{ext}"
    assert user_files.lims_output_dir is None
    assert user_files.gtc_pattern is None
    assert user_files.idat_pattern is None
    assert user_files.ped is None
    assert user_files.map is None
    assert user_files.bed is None
    assert user_files.bim is None
    assert user_files.fam is None


def test_user_files_gtc_pattern():
    # GIVEN-WHEN: we load the user files gtc_pattern with a GTC pattern
    # THEN: it works
    UserFiles(gtc_pattern="{test}.gtc")

    with pytest.raises(ValidationError):
        # GIVEN-WHEN: we load the user files gtc_pattern with a GTC file (no pattern)
        # THEN: it raises a ValidationError
        UserFiles(gtc_pattern="test.gtc")

    with pytest.raises(ValidationError):
        # GIVEN-WHEN: we load the user files gtc_pattern with a TXT patter (not gtc file)
        # THEN: it raises a ValidationError
        UserFiles(gtc_pattern="{test}.txt")


def test_user_files_ped_map():
    # GIVEN-WHEN: we load the user files ped & map with a PED & MAP file
    # THEN: it works
    UserFiles(ped="test.ped", map="test.map")

    with pytest.raises(ValidationError):
        # GIVEN-WHEN: we load the user files ped only
        # THEN: it raises a ValidationError
        UserFiles(ped="test.ped")

    with pytest.raises(ValidationError):
        # GIVEN-WHEN: we load the user files map only
        # THEN: it raises a ValidationError
        UserFiles(map="test.map")


def test_user_files_bed_bim_fam():
    # GIVEN-WHEN: we load the user files bed, bim, & fam
    # THEN: it works
    UserFiles(bed="test.bed", bim="test.bim", fam="test.fam")

    with pytest.raises(ValidationError):
        # GIVEN-WHEN: we load the user files bed & bim only
        # THEN: it raises a ValidationError
        UserFiles(bed="test.bed", bim="test.bim")

    with pytest.raises(ValidationError):
        # GIVEN-WHEN: we load the user files bed & fam only
        # THEN: it raises a ValidationError
        UserFiles(bed="test.bed", fam="test.fam")

    with pytest.raises(ValidationError):
        # GIVEN-WHEN: we load the user files bim & fam only
        # THEN: it raises a ValidationError
        UserFiles(bim="test.bim", fam="test.fam")


def test_idat():
    # GIVEN-WHEN: we load the user files idat_pattern red and green with a idat pattern
    # THEN: it works
    Idat(red="{test}_Red.idat", green="{test}_Grn.idat")

    with pytest.raises(ValidationError):
        # GIVEN-WHEN: we load the user files idat_pattern with no inputs
        # THEN: it raises a ValidationError
        Idat()

    with pytest.raises(ValidationError):
        # GIVEN-WHEN: we load the user files idat_pattern red only
        # THEN: it raises a ValidationError
        Idat(red="{test}_Red.idat")

    with pytest.raises(ValidationError):
        # GIVEN-WHEN: we load the user files idat_pattern red and green but red is not a pattern
        # THEN: it raises a ValidationError
        Idat(red="test_Red.idat", green="{test}_Grn.idat")

    with pytest.raises(ValidationError):
        # GIVEN-WHEN: we load the user files idat_pattern red and green but red is not a idat file
        # THEN: it raises a ValidationError
        Idat(red="{test}_Red.txt", green="{test}_Grn.idat")


def test_software_params_defaults():
    # GIVEN-WHEN: we load the software params with default
    software_params = SoftwareParams()

    # THEN: we get the defaults
    assert software_params.sample_call_rate_1 == 0.8
    assert software_params.snp_call_rate_1 == 0.8
    assert software_params.sample_call_rate_2 == 0.95
    assert software_params.snp_call_rate_2 == 0.95
    assert software_params.ld_prune_r2 == 0.1
    assert software_params.maf_for_ibd == 0.2
    assert software_params.maf_for_hwe == 0.05
    assert software_params.ibd_pi_hat_min == 0.05
    assert software_params.ibd_pi_hat_max == 1.0
    assert software_params.dup_concordance_cutoff == 0.95
    assert software_params.contam_threshold == 0.1
    assert software_params.contam_population == "AF"
    assert software_params.pi_hat_threshold == 0.2
    assert software_params.autosomal_het_threshold == 0.1
    assert software_params.strand == "top"


def test_workflow_params_defaults():
    # GIVEN-WHEN: we load the workflow params with default
    workflow_params = WorkflowParams()

    # THEN: we get the defaults
    assert workflow_params.expected_sex_col_name == "Expected_Sex"
    assert workflow_params.remove_contam is True
    assert workflow_params.remove_sex_discordant is True
    assert workflow_params.remove_rep_discordant is True
    assert workflow_params.remove_unexpected_rep is True
    assert workflow_params.minimum_pop_subjects == 50
    assert workflow_params.control_hwp_threshold == 50


@pytest.mark.parametrize("value", [-0.01, 0, 1.01])
def test_call_rate_param_limits(value):
    """Make sure call rate validation is working.

    Again more a sanity check that pydantic is behaving as expected. Sample
    and SNP call rate filters should between 0-1. Values ≤0 or >1 should
    raise a ``pydantic.ValidationError``.
    """
    # GIVEN-WHEN: we try to set software params outside of the default ranges

    # THEN: they all raise ValidationErrors
    with pytest.raises(ValidationError):
        SoftwareParams(sample_call_rate_1=value)

    with pytest.raises(ValidationError):
        SoftwareParams(snp_call_rate_1=value)

    with pytest.raises(ValidationError):
        SoftwareParams(sample_call_rate_2=value)

    with pytest.raises(ValidationError):
        SoftwareParams(snp_call_rate_2=value)
