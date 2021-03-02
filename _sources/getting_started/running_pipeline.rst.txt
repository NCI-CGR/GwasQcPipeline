Running the Pipeline
====================

There are three phases: configuration, pre-flight checks, and running snakemake or submitting to a cluster.

Configuration
-------------

First, generate a configuration file by running the ``cgr config`` tool. You will be asked for a project name and the location of the sample sheet (manifest). We will then populate the ``config.yml`` file in your current working directory with our default settings

.. code-block:: bash

    $ cgr config  # You will be prompted for the Project Name and Sample Sheet path.

    or

    $ cgr config --project-name "Test Project" --sample-sheet <path to sample sheet>.csv

You should check/edit the ``config.yml`` file. If you are not running on ``CGEMS/CCAD`` you will need to update paths for reference and user files.

.. warning::
    We will raise an error if the sample sheet does not exist.

Pre-flight Checks
-----------------

One of our goals is to fail fast if there are problems. A common source of problems are input reference/user files. We provide the ``cgr pre-flight`` command, which tries to validate all input files to make sure that they are (1) present, (2) readable, and (3) complete.

.. code-block:: bash

    # Everything looks good
    $ cgr pre-flight
    Sample Sheet OK (sample-sheet-name.csv)
    BPM OK (bpm-file-name.bpm)
    VCF OK (vcf-file-name.vcf.gz)
    VCF.TBI OK (tbi-file-name.vcf.gz.tbi)
    Processing GTC files
      [#################################] 100%
    7,231 GTC Files OK.

    or

    # Missing a few GTC files
    $ cgr pre-flight
    Sample Sheet OK (sample-sheet-name.csv)
    BPM OK (bpm-file-name.bpm)
    VCF OK (vcf-file-name.vcf.gz)
    VCF.TBI OK (tbi-file-name.vcf.gz.tbi)
    Processing GTC files
      [#################################] 100%
    There was a problem with these GTC Files:
      FileNotFound:
        - missing-gtc-file1.gtc
        - missing-gtc-file2.gtc

.. warning::
    - These checks can take a while if you have a lot of files.
    - We have not implemented checks for ``idat`` files.


Running Snakemake
-----------------

We use snakemake_ to orchestrate the ``GwasQcPipeline``. We provide a convenience wrapper, ``cgr snakemake``, for interacting/running snakemake locally. This wrapper adds ``-s <path to GwasQcPipeline Snakefile>`` to the snakemake command.

.. _snakemake: https://snakemake.readthedocs.io/en/stable/

To run ``GwasQcPipeline`` locally, you could run:

.. code-block:: bash

    $ cgr snakemake --cores 2 -k --use-conda


Submitting to a cluster
-----------------------

We provide the ``cgr submit`` command to easily submit to a cluster. For NCI users, we support ``CGEMS/CCAD`` and ``Biowulf``. External users will need to create their own `cluster profile`_.

.. _`cluster profile`: https://github.com/snakemake-profiles/doc

.. code-block:: bash

    # Running on CGEMS/CCAD
    $ cgr submit --cgems

    or

    # Running on Biowulf
    $ cgr submit --biowulf

    or

    # Running on external cluster
    $ cgr submit --profile /path/to/my/cluster_profile
    # See https://github.com/snakemake-profiles/doc for suggestions
