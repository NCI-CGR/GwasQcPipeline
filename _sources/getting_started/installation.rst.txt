Installing GwasQcPipeline
=========================

Installing on CGEMS/CCAD
------------------------

These are the recommended installation instructions for ``CGEMS/CCAD``.

Create a ``conda`` environment (python=3.8):

.. code-block:: bash

    $ module load miniconda/4.8.3
    $ conda create -n GwasQcPipeline -y python=3.8 pip
    $ conda activate GwasQcPipeline

Install the current release of the ``GwasQcPipeline``:

.. code-block:: bash

    $ pip install https://storage.googleapis.com/gwasqc/releases/cgr_gwas_qc-1.0.0a1-py3-none-any.whl
    $ cgr --help  # Should provide help information for running the GwasQcPipeline.

After the initial installation, to use ``GwasQcPipeline``:

.. code-block:: bash

    $ module load miniconda/4.8.3
    $ conda activate GwasQcPipeline
    $ cgr --help

To update the ``GwasQcPipeline`` in the future:

.. code-block:: bash

    $ module load miniconda/4.8.3
    $ conda activate GwasQcPipeline
    $ pip install --force-reinstall # <url to latest *.whl>
    # See https://github.com/NCI-CGR/GwasQcPipeline/releases for a list of releases


Installing on other systems
---------------------------

The ``GwasQcPipeline`` requires ``conda`` to run. We suggest you follow these directions to install Miniconda_. Once you have ``conda`` installed, you need to create a ``conda`` environment and install the ``GwasQcPipeline``.

.. _Miniconda: https://conda.io/projects/conda/en/latest/user-guide/install/index.html

.. code-block:: bash

    $ conda create -n GwasQcPipeline -y python=3.8 pip
    $ conda activate GwasQcPipeline
    $ pip install https://storage.googleapis.com/gwasqc/releases/cgr_gwas_qc-1.0.0a1-py3-none-any.whl
    $ cgr --help  # Should provide help information for running the GwasQcPipeline.


To use ``GwasQcPipeline`` first activate your ``conda`` environment:

.. code-block:: bash

    $ conda activate GwasQcPipeline
    $ cgr --help

And to update:

.. code-block:: bash

    $ conda activate GwasQcPipeline
    $ pip install --force-reinstall # <url to latest *.whl>
    # See https://github.com/NCI-CGR/GwasQcPipeline/releases for a list of releases
