Running the Pipeline
====================

Running the pipeline consists of three phases: configuration, pre-flight checks, running snakemake.

Configuration
-------------

To start a new project simply run the ``cgr config`` tool. This will ask you for a project name and location of the sample manifest. It will then populate a ``config.yml`` in your current working directory with common settings for ``cgems``::

    $ cgr config  # if no options are given then you will be prompted for Project Name and Manifest path.

    or

    $ cgr config --project-name "Test Project" --sample-sheet /DCEG/CGF/Infinium/SampleSheets/<file_name>.csv

.. warning::
    The sample sheet must exists or we raise and error.

Now you should have a new file in your working directory called ``config.yml``.

Pre-flight Checks
-----------------

One goal of the pipeline is to fail fast if there are problems. One potential source of issues are user provided files. The ``cgr pre-flight`` command runs a set of data validations on user provided input files to make sure that they are (1) present, (2) readable, and (3) complete.

.. warning::
    This is fully implemented yet.

.. todo::
    Adds docs after implementation of ``pre-flight``.

Running Snakemake
-----------------

Currently, the only way to run the pipeline is locally. We have implemented a thin wrapper around ``snakemake`` called ``cgr snakemake``. You can run ``snakemake`` directly, but you must provide the full path the ``Snakefile`` which can be found by running ``cgr snakemake --help``.

To run the workflow you would do something similar to::

    $ cgr snakemake --cores 2 -k --use-conda

.. warning::
    The current ``Snakefile`` is just a small example to get things working. It does not actually do anything useful.
