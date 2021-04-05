import pytest

from cgr_gwas_qc.models.config import Config
from cgr_gwas_qc.workflow.scripts import plot_autosomal_heterozygosity


@pytest.mark.real_data
def test_plot(real_config: Config, agg_population_qc, tmp_path):
    plot_autosomal_heterozygosity.main(
        agg_population_qc, real_config.software_params.autosomal_het_threshold, tmp_path
    )
