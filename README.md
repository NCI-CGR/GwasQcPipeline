[![read-the-docs](https://img.shields.io/badge/Check%20Out-The%20Docs-blue.svg)](https://nci-cgr.github.io/GwasQcPipeline/)
[![docs-by-sphinx-doc](https://img.shields.io/badge/Docs%20by-Sphinx-1f425f.svg)](https://www.sphinx-doc.org/)

[![default-branch-status](https://github.com/NCI-CGR/GwasQcPipeline/actions/workflows/python-package.yml/badge.svg)](https://github.com/NCI-CGR/GwasQcPipeline/actions/workflows/python-package.yml)
[![built-with-poetry](https://img.shields.io/badge/Built%20with-Poetry-1f425f.svg)](https://python-poetry.org/)

[![snakemake](https://img.shields.io/badge/snakemake-6.4.1-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
[![python](https://img.shields.io/badge/python-3.8-brightgreen.svg)](https://www.python.org/)

[![testing-with-pytest](https://img.shields.io/badge/pytest-enabled-brightgreen.svg)](https://docs.pytest.org/en/stable/)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![linting-with-flake8](https://img.shields.io/badge/flake8-enabled-brightgreen.svg)](https://flake8.pycqa.org/en/latest/)
[![typing-with-mypy](https://img.shields.io/badge/mypy-enabled-brightgreen.svg)](https://mypy.readthedocs.io/en/stable/index.html)


[![code-style-black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![code-style-snakefmt](https://img.shields.io/badge/code%20style-snakefmt-000000.svg)](https://github.com/snakemake/snakefmt)


# GwasQcPipeline

QC pipeline for Illumina SNP Array data generated at [CGR](https://dceg.cancer.gov/about/organization/cgr).

## Usage

If you use this workflow in a paper, please cite the URL of this repository (https://github.com/NCI-CGR/GwasQcPipeline).

The documentation of this pipeline is at https://nci-cgr.github.io/GwasQcPipeline/

### LOG
- add `sex_chr_included` parameter to config.yml. If `false` sex concordance check step is skipped.
- F, M and U are plotted by plot_chrx_inbreeding.py
- Add plink case/control gwas
