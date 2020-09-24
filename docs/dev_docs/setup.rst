Set-up Development Environment
==============================

Package Manager and Build Tool (poetry)
---------------------------------------

This project uses poetry_ as a package manager and build tool. Poetry is a modern python build tool that uses the `pyproject.toml` format to track dependencies and build settings. To install poetry on your system run::

    $ curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python

.. _poetry: https://python-poetry.org/

Once installed you can set up a development environment by running the following::

    $ git clone http://10.133.130.114/jfear/GwasQcPipeline.git
    $ cd GwasQcPipeline
    $ git checkout restructure     # Switch to the active dev branch
    $ poetry install               # Install all dependencies
    $ poetry shell                 # Activate poetry virtual environment

    # Make sure everything is working
    $ make -C docs html                   # This will build documentation
    $ pytest                              # This will run the test suite
    $ mypy                                # This will run type checking
    $ flake8 <path>                       # This will run python linting

    # Build python sdist and wheel
    $ poetry build

This repository also has git pre-commit hooks that need to be activated. Each hook builds a virtual environment, so the first time running will take a little while::

    $ pre-commit install             # Installs the hooks

Now everytime you commit files it will run the required set of tools and make modifications or flag issues found when running these tools.

Workflow Environments
---------------------

The workflow contains a number of compiled 3rd party packages. In order to run the workflow you need to have these tools installed using one of: environment modules, conda, or singularity. For local development conda is the easiest to set-up.

Setting up conda (Miniconda)
::::::::::::::::::::::::::::

We suggest that you use the Miniconda_ distribution which can be installed using the following::

    # Download and install miniconda into your home directory
    $ wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh -O ~/miniconda.sh \
        && /bin/bash ~/miniconda.sh -b  \
        && rm ~/miniconda.sh

    $ conda init --all              # Adds conda to your path
    $ exit                          # You need to restart the terminal for this to take effect
    $ conda update -n base conda    # make conda is working and is up to date

.. _Miniconda: https://docs.conda.io/en/latest/miniconda.html
