Installing GwasQcPipeline
=========================

.. _installation:

Installing on CGEMS/CCAD
------------------------

These are the recommended installation instructions for ``CGEMS/CCAD``.

Create a ``conda`` environment (python=3.8):

.. code-block::

    $ module load miniconda/4.8.3
    $ conda create -n GwasQcPipeline -y python=3.8 pip
    $ conda activate GwasQcPipeline

Install the current release of the ``GwasQcPipeline``:

.. code-block::
    :substitutions:

    $ pip install |pkgurl|
    $ cgr --help  # Should provide help information for running the GwasQcPipeline.

After the initial installation, to use ``GwasQcPipeline``:

.. code-block::

    $ module load miniconda/4.8.3
    $ conda activate GwasQcPipeline
    $ cgr --help

To update the latest version of ``GwasQcPipeline``:

.. code-block::
    :substitutions:

    $ module load miniconda/4.8.3
    $ conda activate GwasQcPipeline
    $ pip install --force-reinstall |pkgurl|
    # See https://github.com/NCI-CGR/GwasQcPipeline/releases for a list of releases

Installing on ccad2
-------------------------

- Install miniconda and then create a GwasQcPipeline production environment
.. code-block::

   $ mkdir /scratch/myfolder/GwasQcPipeline_v1.2
   $ cd /scratch/myfolder/GwasQcPipeline_v1.2
   $ wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
   $ bash Miniconda3-py39_4.12.0-Linux-x86_64.sh -p /scratch/myfolder/GwasQcPipeline_v1.2/conda -b
   $ source conda/bin/activate base
   (base) $ conda update -n base -c defaults conda
   (base) $ conda install -n base -c conda-forge mamba
   (base) $ conda create -n GwasQcPipeline -y python=3.8 pip
   (base) $ conda deactivate

Next install the latest version of the GwasQcPipeline environment:
.. code-block::
   $ source /scratch/myfolder/GwasQcPipeline_v1.2/conda/bin/activate GwasQcPipeline 
   (GwasQcPipeline) $ pip install https://github.com/NCI-CGR/GwasQcPipeline/releases/download/v1.2.0/cgr_gwas_qc-1.2.0-py3-none-any.whl
   (GwasQcPipeline) $ cgr --help
   (GwasQcPipeline) $ cgr version
   v1.2.0

Installing on NIH Biowulf
-------------------------

Running the workflow on Biowulf requires installing miniconda and mamba. A detailed
description is provided on the `Biowulf python env website`_.

.. _`Biowulf python env website`: https://hpc.nih.gov/apps/python.html#envs


.. code-block::

   $ wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.11.0-Linux-x86_64.sh

   # Start an interactive session to access tmp folder. Otherwise installation will fail.
   $ sinteractive --mem=20g --gres=lscratch:20
   $ bash Miniconda3-py39_4.11.0-Linux-x86_64.sh -p /<location of miniconda installation>/conda -b
   $ source /<location of miniconda installation>/conda/bin/activate base
   $ conda update -n base -c defaults conda
   $ conda install -n base -c conda-forge mamba
   $ conda create -n GwasQcPipeline -y python=3.8 pip
   $ conda deactivate

Next install the latest version of the GwasQcPipeline environment:

.. code-block::

   $ source /<location of miniconda installation>/conda/bin/activate GwasQcPipeline
   $ pip install --force-reinstall  https://github.com/NCI-CGR/GwasQcPipeline/releases/download/v1.2.0/cgr_gwas_qc-1.2.0-py3-none-any.whl
   $ cgr --help
   $ cgr version

Installing on other systems
---------------------------

The ``GwasQcPipeline`` requires ``conda`` to run.
We suggest you first install Miniconda_.
Once you have ``conda`` installed, you need to create a ``conda`` environment and install the ``GwasQcPipeline``.

.. code-block::
    :substitutions:

    $ conda create -n GwasQcPipeline -y python=3.8 pip
    $ conda activate GwasQcPipeline
    $ pip install |pkgurl|
    $ cgr --help  # Should provide help information for running the GwasQcPipeline.


To use ``GwasQcPipeline`` first activate your ``conda`` environment:

.. code-block::

    $ conda activate GwasQcPipeline
    $ cgr --help

And to update:

.. code-block::
    :substitutions:

    $ conda activate GwasQcPipeline
    $ pip install --force-reinstall |pkgurl|
    # See https://github.com/NCI-CGR/GwasQcPipeline/releases for a list of releases
