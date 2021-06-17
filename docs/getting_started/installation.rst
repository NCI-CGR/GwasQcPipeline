Installing GwasQcPipeline
=========================

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


Installing on other systems
---------------------------

The ``GwasQcPipeline`` requires ``conda`` to run.
We suggest you follow these directions to install Miniconda_.
Once you have ``conda`` installed, you need to create a ``conda`` environment and install the ``GwasQcPipeline``.

.. _Miniconda: https://conda.io/projects/conda/en/latest/user-guide/install/index.html

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
