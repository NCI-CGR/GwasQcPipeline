from cgr_gwas_qc.cli.config import get_number_samples, get_number_snps


def test_get_number_samples(sample_sheet_file):
    assert 4 == get_number_samples(sample_sheet_file)


def test_get_number_snps(bpm_file):
    assert 2000 == get_number_snps(bpm_file)
