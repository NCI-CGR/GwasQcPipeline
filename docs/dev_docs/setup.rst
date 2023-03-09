.. _dev-setup:

Development Environment Set-up
==============================

Download GwasQcPipeline
-----------------------

First we need to download the development version of GwasQcPipeline.::

    $ git clone --recursive https://github.com/NCI-CGR/GwasQcPipeline.git
    $ cd GwasQcPipeline

.. note::

    Our test data (``tests/data``) is stored in a separate git repository.
    This repository is embedded as a git submodule. The ``--recursive`` tells git to go ahead and download ``tests/data``.


Create a virtual environment (``conda``)
----------------------------------------

We are going to create a ``conda`` virtual environment to store the development environment.
If you need to install ``conda`` see the `Miniconda website`_.

Next you need to setup three channels in your conda config by running the following::

    $ conda config --add channels defaults
    $ conda config --add channels bioconda
    $ conda config --add channels conda-forge

Next, To create our ``conda`` environment run::

    $ conda create -n cgr-dev python=3.8 poetry make -y
    $ conda activate cgr-dev  # This activates the virtual environment

.. _Miniconda website: https://docs.conda.io/en/latest/miniconda.html

.. note::

    ``psutil`` and ``pysam`` are also runtime requirements.
    If you have ``gcc`` installed, they should be automatically built during the installation step below.
    However, if you have problems you may want to install them with ``conda activate cgr-dev && conda install psutil pysam``.

Install dependencies and GwasQcPipeline (``poetry``)
----------------------------------------------------

This project uses poetry_ as a package manager and build tool.
Poetry is a modern python build tool that uses the ``pyproject.toml`` format to track dependencies and build settings.

.. _poetry: https://python-poetry.org/

To install all runtime/development dependencies and GwasQcPipeline itself run::

    $ conda activate cgr-dev      # Make sure we are in our conda environment
    $ poetry install              # Install development and runtime dependencies

Now lets make sure everything is working::

    $ cgr --help                  # This is our main entry point to running the workflow
    $ make -C docs html           # This will build documentation into docs/_build/html
    $ pytest -v                   # This will run the test suite

The main reason we are using poetry is because it makes building python packages easy.
In order to match the new config.yaml version line with the new QwasQcPipeline version
edit and update the ``version = "0.1.2"`` line in the ``pyproject.toml`` file::

    $ poetry build                # Build artifacts are in ./dist

Once the changes are pushed to Github, tag the new version for release. While Github is building
the new release, a drop and drag box will appear for additional assets. Add the new 
cgr_gwas_qc-X.X.X.tar.gz and cgr_gwas_qc-X.X.X-py3-none-any.whl files to the box. 


Install pre-commit hooks for consistent development
---------------------------------------------------

There are a number of tools out there to make coding cleaner and more consistent.
For example, there are code formatters (i.e., ``black``, ``isort``, ``snakefmt``), code linters (i.e., ``flake8``, ``rstcheck``), type checkers (i.e., ``mypy``).
These tools also help catch small mistakes.
This repository has a set of git pre-commit hooks (``.pre-commit-config.yaml``) that will run a suite of tools at each commit.
This helps keeps issues from making it into the code base.

There is a one-time install that you need to setup in your local version of GwasQcPipeline::

    $ pre-commit install             # Installs the hooks
    $ pre-commit run                 # Make sure everything is working

.. note::

    The first time you run pre-commit it needs to download and setup virtual environments for each tool.
    This may take a few minutes.

.. note::

    Tools are only run on files with changes, if this is a fresh clone of the repository then all tools will be skipped.

.. note::

    Now every time you commit files it will run the required set of tools for the staged files.
    If an auto formatter detects a problem it will make the changes, but you will have to re-stage that file.
    This will slow down making commits, but I find the benefits out weight the inconvenience.

.. warning::

    Sometimes pre-commit will keep calling something a problem that you want to ignore.
    For example, ``codespell`` tends to interpret this ``"\nNumber "`` as a spelling error even thought it is really a formatting thing.
    You can skip running al  pre-commit hooks using ``git commit --no-verify``.
    However, make sure it is absolutely necessary!
