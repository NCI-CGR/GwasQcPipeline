#!/usr/bin/env python
# coding: utf-8

"""Test the conversion of BCF to ABF.

For contamination checking we need to have the population level allele
frequencies. We pull these allele frequencies from the 1KG project for each
SNP in the BCF.
"""

import pandas as pd
import pytest

from cgr_gwas_qc.workflow.scripts import vcf2abf


@pytest.mark.regression
def test_vcf2abf_fake_data(fake_data_cache, tmp_path):
    obs_abf = tmp_path / "test.abf.txt"
    vcf2abf.main(
        fake_data_cache / "cgr/bcf/small_genotype.bcf",
        fake_data_cache / "1KG/small_1KG.vcf.gz",
        "AF",
        obs_abf,
    )

    expected_abf = fake_data_cache / "cgr/bcf/small_genotype.abf.txt"

    obs_abf = pd.read_csv(obs_abf, sep="\t")
    expected_abf = pd.read_csv(expected_abf, sep="\t")

    comparision = expected_abf.merge(
        obs_abf, how="left", on="SNP_ID", suffixes=("_expected", "_obs")
    )
    comparision = comparision.dropna()  # excluding the na from the test
    prop_very_different = (abs(comparision.ABF_expected - comparision.ABF_obs) > 1e-6).mean()
    assert 0.01 > prop_very_different  # less than 1% are more than 0.000001 different
