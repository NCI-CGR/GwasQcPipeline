Running the Pipeline
====================

There are three phases to running the |cgr|: configuration, pre-flight checks, and running snakemake or submitting to a cluster.

.. _`cgr-config`:

Creating Configuration
----------------------

Before running the |cgr| we need to create the necessary configuration file (``config.yml``).
We provide a command line utility (``cgr config``) to help generate this file.
To run this utility you need to provide a sample sheet (or CGR LIMs manifest).
Please see :ref:`sample-sheet` to see the sample sheet requirements.
For more details on the available configuration options see :ref:`config-yaml`

.. click:: cgr_gwas_qc.cli.config:typer_click_object
   :prog: cgr config
   :nested: full

.. _`cgr-preflight`:

Pre-Flight File Checks
----------------------

A fundamental design strategy of this project is to fail fast if there is a problem.
A common source of problems are input reference/user files.
We provide the ``cgr pre-flight`` command, which tries to validate all input files to make sure that they are (1) present, (2) readable, and (3) complete.
This command also creates ``cgr_sample_sheet.csv`` which is required by the workflow.
See below for more details about what all ``cgr pre-flight`` does.

Here is example output where checks pass::

    Sample Sheet OK (sample-sheet-name.csv)
    BPM OK (bpm-file-name.bpm)
    VCF OK (vcf-file-name.vcf.gz)
    VCF.TBI OK (tbi-file-name.vcf.gz.tbi)
    Processing GTC files
      [#################################] 100%
    7,231 GTC Files OK.

Here is example output with two missing GTC files::

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

.. attention::
   If the ``config.yml``, sample sheet, or reference files have any issues you must fix them before continuing.
   If you are missing a few sample's IDAT or GTC files you may decide to continue;
   the workflow will automatically exclude these samples.

.. click:: cgr_gwas_qc.cli.pre_flight:typer_click_object
   :prog: cgr pre-flight
   :nested: full

.. _`cgr-snakemake`:

Running the Workflow Locally
----------------------------

We use snakemake_ to orchestrate the |cgr|.
We provide ``cgr snakemake`` to simplify running snakemake_ directly.

.. click:: cgr_gwas_qc.cli.snakemake:typer_click_object
   :prog: cgr snakemake
   :nested: full

.. _`cgr-submit`:

Running the Workflow on a Cluster
---------------------------------

We provide the ``cgr submit`` command to easily submit to a cluster.
We take advantage of snakemake_'s cluster profile system to run on different cluster environments.
For CGR users, we include cluster profiles for ``CGEMS/CCAD``, ``ccad2``, and ``Biowulf``.
For external users will need to create your own `snakemake cluster profile`_.

.. _`snakemake cluster profile`: https://github.com/snakemake-profiles/doc

.. click:: cgr_gwas_qc.cli.submit:typer_click_object
   :prog: cgr submit
   :nested: full

.. note::
   External users on a SLURM or SGE cluster, may just want to modify one of our `profiles <cluster_profiles_>`_.

.. attention::

   **Biowulf users**.
   You may need to adjust ``--time-hr``, ``--local_mem_mb``, and ``--local_tasks``
   if your main job is getting killed by the cluster because of resource limits.

The submission script will create a log file ``./gwas_qc_log.$JOB_ID`` that will have the status of your cluster submission.
Logs for each submitted job can be found in ``./logs/``.
