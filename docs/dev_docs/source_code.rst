.. _source:

Working with Source Code
========================

The |cgr| source code is located in ``src/cgr_gwas_qc``.
This includes the GwasQcPipeline workflow, the ``cgr`` command line utility, and various helper functions.

.. code-block:: bash

   .
   ├── src
   │  └── cgr_gwas_qc
   │     ├── __init__.py
   │     ├── __main__.py          # Points to the ``cgr`` entry point
   │     ├── cli/                 # All code for the ``cgr`` command line interface
   │     ├── cluster_profiles/    # Cluster profiles for CGEMs and Biowuld
   │     ├── conda.py             # Helpers to activate conda environments in a script
   │     ├── config.py            # Configuration manager helper code ``cfg``
   │     ├── exceptions.py        # Custom exceptions for nicer error handling
   │     ├── models/              # Data models for ``config.yml``
   │     ├── parsers/             # Various tools for parsing 3rd party files
   │     ├── paths.py             # Helpers for working with paths
   │     ├── reporting/           # Code for writing the QC report ``.docx``
   │     ├── scripts/             # Legacy compare script
   │     ├── testing/             # Helpers for testing
   │     ├── typing.py            # Custom PathLike data type
   │     ├── validators/          # Code to validate different file formats
   │     ├── version.py           # Gets the workflow version
   │     ├── workflow/            # The GwasQcPipeline workflow
   │     └── yaml.py              # Helpers for working with YAML files


I will briefly discussion a few of the more important or complicated bits.

Configuration Management
------------------------

.. code-block:: bash

   src
   └── cgr_gwas_qc
      ├── ...
      ├── config.py                   # The ``cgr_gwas_qc.config.ConfigMgr``
      ├── models/                     # Pydantic data models defining ``config.yml``
      │  ├── __init__.py
      │  └── config/
      │     ├── __init__.py           # Top Level ``config.yml``
      │     ├── reference_files.py    # Reference Files namespace
      │     ├── software_params.py    # Software Parameters namespace
      │     ├── user_files.py         # User Files namespace
      │     └── workflow_params.py    # Workflow Parameters
      ├── parsers/
      │  ├── sample_sheet.py          # LIMs Manifest and sample sheet parsers
      ├── ...

Configuration management is the most complicated aspect of the entire project.
I wanted to create an easy to use system for all things configuration.
The core of this system is ``cgr_gwas_qc.config.ConfigMgr``.
You will see this referred to as ``cfg`` throughout the workflow, tests, and some scripts.
``cfg`` provides easy access to the values stored in the ``config.yml` as a dictionary at ``cfg.config``, the ``cgr_sample_sheet.csv`` as a data frame at ``cfg.ss``, and to a number of helpers.

Essentially, I am trying to solve two problems with ``cfg``.

1. **Resolving file paths.**
snakemake_ makes many assumptions about the locations of inputs/outputs, snakefiles, conda environments, and scripts.
Traditionally this works as expected when you have your workflow in the current working directory.
However, GwasQcPipeline lives in the install location so we need to adjust adjust the paths for snakemake_.
``cgf`` has several helpers to do just that.

.. autoclass:: cgr_gwas_qc.config.ConfigMgr
   :members: SRC_DIR, WORKFLOW_DIR, CONDA_DIR, MODULE_DIR, SCRIPTS_DIR, SUBWORKFLOW_DIR, SNAKEFILE, TEMPLATE_DIR, conda, scripts, subworkflow, docx_template
   :noindex:

2. **Using config/sample sheet values.**
A common pattern in snakemake_ workflows is to fill in wildcards of a file name pattern.
We often want to use values stored in ``cgr_sample_sheet.csv`` to do this.
The ``cgf`` object gives us access to both the ``config.yml`` and ``cgr_sample_sheet.csv``.
And it provides an enhanced version of snakemake_'s ``expand`` function to minimize the amount of code you have to type.

.. autoclass:: cgr_gwas_qc.config.ConfigMgr
   :members: config, ss, cluster_groups, expand
   :noindex:


As discussed in :ref:`config-yaml`, the other part of the configuration system are the data models that describe ``config.yml``.
These models are located in ``cgr_gwas_qc.models.config``.
Here I am using pydantic_ to describe the ``config.yml``.
Essentially, that means I create a very basic class to define each namespace.
I was not sure if I was going to have other data models, so that is why these live in the ``models/config`` subfolder.
I chose the ``models`` terminology because that is commonly used in the web application space.

To add a new config option, you need to do the following:

1. Decide which name space you want the value to reside.
2. Add that value to the correct class in ``models/config``.
3. You can give the option a default value in the pydantic class (suggested for ``software_params`` and ``workflow_params``), or you can add a default value to ``cli/config.py`` (suggested for ``reference_files`` and ``user_files``).
4. You may need to adjust tests if you require a value but do not provide a default.

GwasQcPipeline Workflow
-----------------------

.. code-block:: bash

   src
   └── cgr_gwas_qc
      ├── workflow
      │  ├── __init__.py
      │  ├── conda/                       # Definition of conda environments used by the workflow
      │  │  ├── eigensoft.yml
      │  │  ├── graf.yml
      │  │  ├── illuminaio.yml
      │  │  ├── king.yml
      │  │  ├── pandoc.yml
      │  │  ├── plink2.yml
      │  │  └── verifyidintensity.yml
      │  ├── modules/                     # Rules for different software packages
      │  │  ├── eigensoft.smk
      │  │  ├── graf.smk
      │  │  ├── plink.smk
      │  │  └── thousand_genomes.smk
      │  ├── scripts/                     # Custom python scripts used by the workflow
      │  │  ├── __init__.py
      │  │  └── ...
      │  ├── Snakefile                    # The main workflow, this is what is called when a sub-workflow is not specified
      │  └── sub_workflows/               # The sub-workflows
      │     ├── contamination.smk
      │     ├── delivery.smk
      │     ├── entry_points.smk
      │     ├── sample_qc.smk
      │     └── subject_qc.smk

The snakemake_ workflow lives in ``cgr_gwas_qc.workflow``.
When you run ``cgr snakemake`` or ``cgr submit`` without specifying a subworkflow you run ``workflow/Sankefile``.
If you use the ``--subworkflow <name>`` then you will run one of ``workflow/subworkflow``

Besides sub-workflow, I also have a number of software modules located in `workflow/modules`.
My goal here was to group all of the rules for each 3rd party software to make maintenance easier.
I also tried to follow the key programming principle, only do something once!

.. note::
   There are a few exceptions to all 3rd party tools are in a module.

   - I had to embed a PLINK command in ``workflow/scripts/plink_merge.py``.
   - I had to embed verifyidintensity in ``workflow/scripts/verifyidintensity.py`
   - I never created a King module, so all that code is in ``workflow/sub_workflows/sample_qc.smk``

Almost all custom code are in ``cgr_gwas_qc.workflow.scripts``.
These scripts are designed to be called directly from snakemake using the ``scripts:`` directive.
Some of them can also be run directly from the command line.
I tried to use descriptive names, but in general it is probably easier to go to the corresponding sub-workflow and see which script is referenced.


QC Report Generation
--------------------

.. code-block:: bash

   src
   └── cgr_gwas_qc
      ├── reporting
      │  ├── __init__.py
      │  ├── constants.py
      │  ├── qc_exclusions.py
      │  ├── sample_qc.py
      │  ├── subject_qc.py
      │  ├── templates
      │  │  ├── __init__.py
      │  │  ├── cgr_reference.docx
      │  │  ├── qc_report.md.j2
      │  │  └── summary_stats.txt.j2
      │  └── templating.py

One major change from the legacy workflow is that we moved the QC report to a markdown file that is converted to ``docx`` at runtime.
The text of the QC report is in ``cgr_gwas_qc/reporting/templates/qc_report.md.j2``.
This file is in the jinja2_ format.
jinja2_ is a template framework.
It is mostly used to dynamically build HTML pages, but here we use it to build our markdown report.
It allows for value substitutions, flow control (if/else), and loops.
To build the QC report we create a giant dictionary with all of the values and then pass this "payload" to jinja2_.

The report is pretty complicated and has a lot of logic to generate all of the necessary counts.
The workflow script ``workflow/scripts/qc_report.py`` uses helpers in the reporting folder to try to keep things clean.
There are three major sections: Sample QC, Subject QC, QC Exclusions (a.k.a Table 4).
The values in these three sections are calculated in ``reporting/sample_qc.py``, ``reporting/subject_qc.py``, and ``reporting/qc_exclusions.py`` respectively.


Command Line Interface (CLI)
----------------------------

.. code-block:: bash

   src
   └── cgr_gwas_qc
      ├── __main__.py
      ├── cli
      │  ├── __init__.py
      │  ├── cgems_copy.py
      │  ├── config.py
      │  ├── pre_flight.py
      │  ├── snakemake.py
      │  ├── submit.py
      │  └── version.py

The ``cgr`` command line interface lives in ``cgr_gwas_qc/cli``.
Each of the sub-commands are in their own python file.
I am using the python CLI framework called typer_.
It is an easy way to build up nice command line tools.
To add a new command you would create a new python script in ``cgr_gwas_qc/cli`` and register the command with ``cgr`` by adding it to ``cli/__init__.py``


NIH Cluster Profiles
--------------------

.. code-block:: bash

   src
   └── cgr_gwas_qc
      ├── cluster_profiles
      │  ├── __init__.py
      │  ├── biowulf
      │  │  ├── __init__.py
      │  │  ├── config.yaml
      │  │  ├── jobscript.sh
      │  │  ├── status.py
      │  │  └── submit.py
      │  ├── cgems
      │  │  ├── __init__.py
      │  │  ├── cgems_jobscript.sh
      │  │  ├── cgems_status.py
      │  │  ├── cgems_submit.py
      │  │  └── config.yaml
      │  ├── cluster.yaml
      │  └── snakemake.sh

snakemake_ provides a plugin architecture to tell it how to interact with your cluster environment (`their docs <https://github.com/snakemake-profiles/doc>`__).
This consists of 3 major parts: a job-script template, a submission script, and a status check script.
For each rule that is submitted to the cluster, snakemake creates a job.
The submission script is given this job, it can then pull out job's metadata and build a submission script using the job-script template.
snakemake_ will run the status check script repeatedly for each submitted job.
The status script then tells snakemake if the job is "running", "failed", or "success".

I created two cluster profiles, one for Biowulf and one for CGEMs/CCAD.
Code that is shared by both profiles lives in ``cluster_profiles/__init__.py``

- job-script template: ``cluster_profiles/biowulf/jobscript.sh`` and ``cluster_profiles/cgems/cgems_jobscript.sh``
- submit-script template: ``cluster_profiles/biowulf/submit.py`` and ``cluster_profiles/cgems/cgems_submit.py``
- status-script template: ``cluster_profiles/biowulf/status.py`` and ``cluster_profiles/cgems/cgems_status.py``

When you run ``cgr submit --cgems`` or ``cgr submit --biowulf`` we create a submission script for the main snakemake task that is saved in ``.snakemake/GwasQcPipeline_submission.sh``
To create this script we use the jinja2_ template ``cluster_profiles/snakemake.sh``.


Parsers and Validators
----------------------

.. code-block:: bash

   src
   └── cgr_gwas_qc
      ├── parsers
      │  ├── __init__.py
      │  ├── bim.py
      │  ├── bpm.py
      │  ├── common.py
      │  ├── eigensoft.py
      │  ├── graf.py
      │  ├── illumina
      │  │  ├── __init__.py
      │  │  ├── adpc.py
      │  │  └── IlluminaBeadArrayFiles.py
      │  ├── king.py
      │  ├── plink.py
      │  ├── sample_sheet.py
      │  ├── vcf.py
      │  └── verifyidintensity.py
      ├── validators
      │  ├── __init__.py
      │  ├── bgzip.py
      │  ├── bpm.py
      │  ├── gtc.py
      │  ├── idat.py
      │  └── sample_sheet.py

This workflow uses a lot of different file types.
Many of them are binary while others are plain text in a various formats.
I created a series of parsers and validators for many of the different file types.
I centralized them here again to only do something once!

The file ``cgr_gwas_qc/parsers/illumina/IlluminaBeadArrayFiles.py`` was downloaded from Illumina's github page.
All of the other parsers/validators I wrote.
See :ref:`API Reference` for more details.
