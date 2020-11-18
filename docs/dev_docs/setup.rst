Setting-up the Development Environment
======================================

The GwasQcPipeline is a python based workflow that uses snakemake for orchestration. There are a number of ways to setup the development environment, but bellow is our recommend method.

Download GwasQcPipeline
-----------------------

First we need to download the development version of GwasQcPipeline. Major development is occurring on the ``restructure`` branch::

    $ git clone -b restructure --recursive http://10.133.130.114/jfear/GwasQcPipeline.git
    $ cd GwasQcPipeline

.. note::

    Our test data (``tests/data``) is stored in a separate git repository. This repository is embedded as a git submodule. The ``--recursive`` tells git to go ahead and download ``tests/data``.


Create a virtual environment (``conda``)
----------------------------------------

We recommend you to use virtual environments. There are a number of ways to do this, but since ``conda`` is needed to run the workflow we suggest you just use a ``conda`` environment. If you need to install ``conda`` see the `Miniconda website`_.

Next you need to add two repositories to your conda config by running the following::

    $ conda config --add channels defaults
    $ conda config --add channels bioconda
    $ conda config --add channels conda-forge

Next, To create our ``conda`` environment run::

    $ conda create -n GwasQcPipeline python=3.8 poetry psutil pysam make -y
    $ conda activate GwasQcPipeline        # Activates virtual environment

.. _Miniconda website: https://docs.conda.io/en/latest/miniconda.html

.. note::

    ``psutil`` and ``pysam`` are runtime requirements. You can install both using poetry/pip (below), but then ``psutil`` requires ``gcc`` to be installed. Similarly, ``pysam`` is sometimes problematic so it is best just to use the ``conda`` version.

Install dependencies and GwasQcPipeline (``poetry``)
----------------------------------------------------

This project uses poetry_ as a package manager and build tool. Poetry is a modern python build tool that uses the ``pyproject.toml`` format to track dependencies and build settings.

.. _poetry: https://python-poetry.org/

To install all runtime/development dependencies and GwasQcPipeline itself run::

    $ conda activate GwasQcPipeline    # Make sure we are in our conda environment
    $ poetry install                   # Install development and runtime dependencies

Now lets make sure everything is working::

    $ make -C docs html                   # This will build documentation into docs/_build/html
    $ pytest -v                           # This will run the test suite
    $ mypy                                # This will run type checking
    $ gtc2plink --help                    # Is one of this projects command line tools

The main reason we are using poetry is because it makes building python packages easy. This won't be useful until we are ready to deploy. However, if you wanted to build a python package that you could push to pypi or email someone you would run::

    $ poetry build                        # Build artifacts are in ./dist

Install pre-commit hooks for consistent development
---------------------------------------------------

There are a number of tools out there to make coding cleaner and more consistent. For example, there are code formatters (i.e., ``black``, ``isort``, ``snakefmt``), code linters (i.e., ``flake8``, ``rstcheck``), type checkers (i.e., ``mypy``). These tools also help catch small mistakes. This repository has a set of git pre-commit hooks (``.pre-commit-config.yaml``) that will run a suite of tools at each commit. This helps keeps issues from making it into the code base.

There is a one-time install that you need to setup in your local version of GwasQcPipeline::

    $ pre-commit install             # Installs the hooks
    $ pre-commit run                 # Make sure everything is working

.. note::

    The first time you run pre-commit it needs to download and setup virtual environments for each tool. This may take a few minutes.

.. note::

    Tools are only run on files with changes, if this is a fresh clone of the repository then all tools will be skipped.

.. warning::

    Now every time you commit files it will run the required set of tools for the staged files. If an auto formatter detects a problem it will make the changes, but you will have to re-stage that file.

    This will slow down making commits, but I find the benefits out weight the inconvenience. If you want to force something and ignore the pre-commit hooks you can always run ``git commit --no-verify``. However, make sure it is absolutely necessary!
