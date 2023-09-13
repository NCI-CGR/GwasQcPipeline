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

### Deploying with Docker
Build docker image from within GwasQcPipeline directory
```
docker build -t gwas_qc_pipe .
```

Test docker image if you have test data
```
docker run -v $(pwd):/home/data -i -t gwas_qc_pipe snakemake -k --use-conda -npr
```

# Installing GwasQcPipeline on ccad2
- Install miniconda and then create GwasQcPipeline production environment
```
$ mkdir /scratch/myfolder/GwasQcPipeline_v1.2
$ cd /scratch/myfolder/GwasQcPipeline_v1.2
$ wget https://repo.anaconda.com/miniconda/ Miniconda3-py39_4.12.0-Linux-x86_64.sh
$ bash Miniconda3-py39_4.12.0-Linux-x86_64.sh -p /scratch/myfolder/GwasQcPipeline_v1.2/conda -b
$ source conda/bin/activate base
(base) $ conda update -n base -c defaults conda
(base) $ conda install -n base -c conda-forge mamba
(base) $ conda create -n GwasQcPipeline -y python=3.8 pip
(base) $ conda deactivate
```
- Install GwasQcPipeline source code
```
$ source /scratch/myfolder/GwasQcPipeline_v1.2/conda/bin/activate GwasQcPipeline 
(GwasQcPipeline) $ pip install https://github.com/NCI-CGR/GwasQcPipeline/releases/download/v1.2.0/cgr_gwas_qc-1.2.0-py3-none-any.whl
(GwasQcPipeline) $ cgr --help
(GwasQcPipeline) $ cgr version
v1.2.0
```
## Usage
```cgr submit --ccad2 --local-mem-mb 8000```

### LOG v1.2.0
- add ccad2-slurm profile

