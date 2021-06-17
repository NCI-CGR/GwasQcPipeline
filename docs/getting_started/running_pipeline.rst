Running the Pipeline
====================

There are three phases: configuration, pre-flight checks, and running snakemake or submitting to a cluster.

Configuration
-------------

First, generate a configuration file by running the ``cgr config`` tool.
You must provide the location of the sample sheet (or LIMs manifest).
We will then populate the ``config.yml`` file with our default settings.

.. code-block::

    $ cgr config  # You will be prompted for the Sample Sheet path.

    or

    $ cgr config --sample-sheet <path to sample sheet>.csv

You will then need to edit/check the ``config.yml`` file to make sure paths and settings are correct.

.. warning::
   We will raise an error if the sample sheet does not exist.

.. note::
   See ``cgr config --help`` for more options.

CGEMS/CCAD
^^^^^^^^^^

For ``CGEMS/CCAD`` users we have a few additional options to streamline config generation.

.. code-block::

  $ cgr config --cgems --sample-sheet <path to sample sheet>.csv

  or

  $ cgr config --cgems-dev --sample-sheet <path to sample sheet>.csv

Both of these options populate ``reference_files`` and ``user_files`` with the standard paths found on ``CGEMS/CCAD``.
The option ``--cgems`` creates the standard production run folder structure and saves the ``config.yml`` there.
The option ``--cgems-dev`` does not create the production run folders, and saves the ``config.yml`` in the current working directory.

At CGR it is also common practice to re-run a project with slightly different settings.
We include a helper command ``cgr cgems-copy`` that will clone the ``config.yml`` and the ``cgr_sample_sheet.csv`` (see Pre-flight Checks) to a new production run folder.

.. code-block::

   $ cd <...>/GSA_Lab_QC/<project name>/builds
   $ ls
   QC_v1_01012021
   QC_v2_01042021

   $ cgr cgems-copy .
   $ ls
   QC_v1_01012021
   QC_v2_01042021
   QC_v3_<today's date>

   $ ls QC_v3_*
   cgr_sample_sheet.csv
   config.yml

Pre-flight Checks
-----------------

A fundamental design strategy of this project is to fail fast if there is a problem.
A common source of problems are input reference/user files.
We provide the ``cgr pre-flight`` command, which tries to validate all input files to make sure that they are (1) present, (2) readable, and (3) complete.

.. code-block::

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
   - These checks can take a while if you have a lot of files, but we do provide ``cgr pre-flight --threads <num threads>`` to allow parallel processing.

.. note::
   We also use pre-flight checks to normalize the sample sheet information used while running the workflow.
   At the CGR our samples sheets are actually manifest files from our LIMs system (INI format).
   So regardless of if you provide a basic sample sheet or a manifest file we generate a normalized sample sheet ``cgr_sample_sheet.csv`` that is used by the workflow.

.. note::
   See ``cgr pre-flight --help`` for more options.

Running Snakemake
-----------------

We use snakemake_ to orchestrate the ``GwasQcPipeline``.
We provide a convenience wrapper, ``cgr snakemake``, for interacting/running snakemake locally.
This wrapper adds ``-s <path to GwasQcPipeline Snakefile>`` to the snakemake command.

.. _snakemake: https://snakemake.readthedocs.io/en/stable/

To run ``GwasQcPipeline`` locally, you could run:

.. code-block::

    $ cgr snakemake --cores 2 -k --use-conda

.. note::
   See ``cgr snakemake --help`` for more options.

Submitting to a cluster
-----------------------

We provide the ``cgr submit`` command to easily submit to a cluster.
We take advantage of snakemake_'s cluster profile system to run on different cluster environments.
For CGR users, we include cluster profiles for ``CGEMS/CCAD`` and ``Biowulf``.
For external users will need to create your own `snakemake cluster profile`_.

.. _`snakemake cluster profile`: https://github.com/snakemake-profiles/doc

.. code-block::

    # Running on an external cluster
    $ cgr submit --profile /path/to/my/cluster_profile --queue <queue name> --submission-cmd <name submission command>

.. note::
   External users on a SLURM or SGE cluster, may just want to modify one of our `profiles <https://github.com/NCI-CGR/GwasQcPipeline/tree/default/src/cgr_gwas_qc/cluster_profiles>`_.

.. note::
   See ``cgr submit --help`` for more options.

CGEMS/CCAD and Biowulf
^^^^^^^^^^^^^^^^^^^^^^

To submit the workflow on CGEMS/CCAD or Biowulf use the following options.

.. code-block::

    # Running on CGEMS/CCAD
    $ cgr submit --cgems

    or

    # Running on Biowulf
    $ cgr submit --biowulf

This will generate a submission script in ``.snakemake/GwasQcPipeline_submission.sh`` and submit this job using ``qsub`` or ``sbatch``.
It is possible to directly edit this file and resubmit without using `cgr submit` again.

.. warning::

   Attention Biowulf users.
   A number of steps are run "locally" as part of the main snakemake job.
   We try to provide sane defaults for jobs, but if your main job is killed by the cluster because of resource limits you may need to adjust ``--time-hr``, ``--local-mem-mb``, and ``--local-tasks``.
