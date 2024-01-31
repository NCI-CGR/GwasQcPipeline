.. _dev-setup:

Development Environment Set-up
==============================


Install ``poetry``
------------------

First we need to download poetry by using the official `Poetry installer`_.::

    $ curl -sSL https://install.python-poetry.org | python -
    $ poetry --version

.. _Poetry installer: https://python-poetry.org/docs/#installation
.. _instructions: https://www.baeldung.com/linux/default-python3

.. note::

    If poetry version is not accessible, check your PATH and ensure poetry's install location is findable.
    If default python version is not 3.8, then follow these `instructions`_.

Create a virtual environment (``conda``)
----------------------------------------

We are going to create a ``conda`` virtual environment to store the development environment.
If you need to install ``conda`` see the `Miniconda website`_. We recomend installing conda locally::

    $ wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh
    $ bash Miniconda3-py38_4.12.0-Linux-x86_64.sh -p miniconda3 -b
    $ echo "export PATH=$(pwd)/miniconda3/bin:$PATH" >> ~/.bashrc
    $ source ~/.bashrc
    $ conda init

Restart you terminal.

Next you need to setup three channels in your conda config by running the following::

    $ export PYTHON_KEYRING_BACKEND=keyring.backends.null.Keyring
    $ conda update -n base -c defaults conda
    $ conda install -n base -c conda-forge mamba
    $ conda config --add channels defaults
    $ conda config --add channels bioconda
    $ conda config --add channels conda-forge

Next, To create our ``conda`` environment run::

    $ conda create -n cgr-dev python=3.8 -y
    $ conda activate cgr-dev  # This activates the virtual environment

.. _Miniconda website: https://docs.conda.io/en/latest/miniconda.html

.. note::

    ``psutil`` and ``pysam`` are also runtime requirements.
    If you have ``gcc`` installed, they should be automatically built during the installation step below.
    However, if you have problems you may want to install them with ``conda activate cgr-dev && conda install psutil pysam``.

Download GwasQcPipeline
-----------------------

Download the development version of GwasQcPipeline.::

    $ git clone --recursive https://github.com/NCI-CGR/GwasQcPipeline.git
    $ cd GwasQcPipeline

.. note::

    Our test data (``tests/data``) is stored in a separate git repository.
    This repository is embedded as a git submodule. The ``--recursive`` tells git to go ahead and download ``tests/data``.


Install dependencies and GwasQcPipeline (``poetry``)
----------------------------------------------------

This project uses poetry_ as a package manager and build tool.
Poetry is a modern python build tool that uses the ``pyproject.toml`` format to track dependencies and build settings.

.. _poetry: https://python-poetry.org/

To install all runtime/development dependencies and GwasQcPipeline itself run::

    $ conda activate cgr-dev      # Make sure we are in our conda environment
    $ poetry env use /PATH/TO/miniconda3/envs/cgr-dev/bin/python # Enable poetry to manage your conda environment
    $ poetry config virtualenvs.path /PATH/TO/miniconda3/envs/cgr-dev #This needs to be full path
    $ poetry env info # This should show that both system and virtual env python is 3.8 and that the venv is conda
    $ poetry install              # Install development and runtime dependencies
    $ cgr version

Now lets make sure everything is working::

    $ cgr --help                  # This is our main entry point to running the workflow
    $ make -C docs html           # This will build documentation into docs/_build/html
    $ poetry run pytest -vvv      # This will run the test suite

The main reason we are using poetry is because it makes building python packages easy.
In order to match the new config.yaml version line with the new QwasQcPipeline version
edit and update the ``version = "1.2.0"`` line in the ``pyproject.toml`` file::

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
This helps keep issues from making it into the code base.

There is a one-time install that you need to setup in your local version of GwasQcPipeline::

    $ pre-commit install             # Installs the hooks
    $ pre-commit run                 # Make sure everything is working

.. note::

    The first time you run pre-commit it needs to download and setup virtual environments for each tool.
    This may take a few minutes.

.. note::

    Tools are only run on files with changes, if this is a fresh clone of the repository then all tools will be skipped.

.. note::

    Now, every time you commit files, it will run the required set of tools for the staged files.
    If an auto formatter detects a problem, it will make the changes, but you will have to re-stage that file.
    This will slow down making commits, but I find the benefits out weight the inconvenience.

.. warning::

    Sometimes pre-commit will keep calling something a problem that you want to ignore.
    For example, ``codespell`` tends to interpret this ``"\nNumber "`` as a spelling error even thought it is really a formatting thing.
    You can skip running all  pre-commit hooks using ``git commit --no-verify``.
    However, make sure it is absolutely necessary!
