import pytest

from cgr_gwas_qc.config import ConfigMgr
from cgr_gwas_qc.workflow.scripts import plot_autosomal_heterozygosity


@pytest.mark.real_data
def test_plot(real_cfg: ConfigMgr, agg_population_qc, tmp_path):
    plot_autosomal_heterozygosity.main(
        agg_population_qc, real_cfg.config.software_params.autosomal_het_threshold, tmp_path
    )
