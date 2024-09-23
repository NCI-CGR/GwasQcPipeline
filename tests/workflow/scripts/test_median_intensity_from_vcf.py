#!/usr/bin/env python
# coding: utf-8

"""Test calculation of median intensity from BCF

We need to assess whether the intensity meets the miniumum threshold. For this purpose, we typically check median IDAT intensity. Here we attempt to test the same with BCF as input file.
"""

import pandas as pd
import pytest
from numpy import isclose

from cgr_gwas_qc.workflow.scripts import median_intensity_from_vcf


@pytest.mark.regression
def test_median_intensity_from_vcf_fake_data(fake_data_cache, tmp_path):
    obs_median_int = tmp_path / "test_median_int.csv"
    median_intensity_from_vcf.main(
        fake_data_cache / "cgr/bcf/small_genotype.bcf", "small_genotype", obs_median_int
    )

    expected_median_int = fake_data_cache / "cgr/bcf/small_genotype_median_int.csv"

    expected_median_int = pd.read_csv(expected_median_int)
    obs_median_int = pd.read_csv(obs_median_int)
    assert isclose(expected_median_int.median_intensity, obs_median_int.median_intensity)
