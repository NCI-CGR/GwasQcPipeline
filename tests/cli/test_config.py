from cgr_gwas_qc.cli.config import get_number_snps


def test_get_number_snps(bpm_file):
    assert 2000 == get_number_snps(bpm_file)
