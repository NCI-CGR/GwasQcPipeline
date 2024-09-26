#!/usr/bin/env python
# coding: utf-8

"""Test the extraction of sample scores from BCF and writing to adpc.bin.

For contamination checking we need to have the sample scores written into adpc.bin to serve as input for verifyIDintensity.
"""

import pytest
from numpy import isclose

from cgr_gwas_qc.parsers.illumina import AdpcReader
from cgr_gwas_qc.workflow.scripts import vcf2adpc


@pytest.mark.regression
def test_vcf2adpc_fake_data(fake_data_cache, tmp_path):
    obs_adpc = tmp_path / "test.adpc.bin"
    vcf2adpc.main(fake_data_cache / "cgr/bcf/small_genotype.bcf", "small_genotype", obs_adpc)

    expected_adpc = fake_data_cache / "cgr/bcf/small_genotype.adpc.bin"

    with AdpcReader(expected_adpc) as expected, AdpcReader(obs_adpc) as obs:
        for exp_row, obs_row in zip(expected, obs):
            # Integers should be the exactly the same
            assert exp_row.x_raw == obs_row.x_raw
            assert exp_row.y_raw == obs_row.y_raw
            assert exp_row.genotype == obs_row.genotype

            # Floats should be really close
            assert isclose(exp_row.x_norm, obs_row.x_norm)
            assert isclose(exp_row.y_norm, obs_row.y_norm)
            assert isclose(exp_row.genotype_score, obs_row.genotype_score)
