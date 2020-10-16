Installing GwasQcPipeline
=========================

Installing on CGEMs
-------------------

These are the recommended installation instructions for ``cgems`` HPC.

Create a ``conda`` environment::

    $ module load miniconda/4.8.3
    $ conda create -n GwasQcPipeline -y python=3.8 pip
    $ conda activate GwasQcPipeline

.. note::
    The GwasQcPipeline requires ``python>=3.8``. On ``cgems`` this is easiest to achieve by using a conda environment.

Install ``GwasQcPipeline`` from the current version on ``cgems``::

    $ pip install /DCEG/CGF/Bioinformatics/Production/fearjm/cgr_gwas_qc-0.1.0-py3-none-any.whl
    $ cgr --help  # Should provide help information for running the GwasQcPipeline.

To use ``GwasQcPipeline``::

    $ module lod miniconda/4.8.3
    $ conda activate GwasQcPipeline
    $ cgr --help

To update the ``GwasQcPipeline`` in the future::

    $ module lod miniconda/4.8.3
    $ conda activate GwasQcPipeline
    $ pip install --force-reinstall /DCEG/CGF/Bioinformatics/Production/fearjm/cgr_gwas_qc-0.1.0-py3-none-any.whl


Installing on other systems
---------------------------

Currently, the install able wheel is only available on ``cgems``. You could copy ``/DCEG/CGF/Bioinformatics/Production/fearjm/cgr_gwas_qc-0.1.0-py3-none-any.whl`` else where and follow the procedures outlined above.

You can download and install from the GitLab repository with the following directions.

Again setup a ``conda`` environment but this time install ``poetry`` instead of ``pip``::

    $ conda create -n GwasQcPipeline -y python=3.8 poetry
    $ conda activate GwasQcPipeline


Clone the development version of GwasQcPipeline. Major development is occurring on the ``restructure`` branch::

    $ git clone -b restructure --recursive http://10.133.130.114/jfear/GwasQcPipeline.git
    $ cd GwasQcPipeline

Install the package using ``poetry``::

    $ poetry install --no-dev

.. note::
    For more details see development docs.
