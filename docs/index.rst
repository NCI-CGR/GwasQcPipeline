GwasQcPipeline Documentation
============================

The |cgr| generates a number of sample level QC metrics followed by several filtering criteria.
Samples are then split into ancestral populations and examined for relatedness, population structure, and genotyping errors.
The full workflow is pictured in :numref:`fig-full-pipeline`.

.. figure:: static/GwasQcPipelineWorkflow.*
   :name: fig-full-pipeline

   The |cgr|.

.. toctree::
   :caption: Getting Started
   :name: getting-started
   :maxdepth: 1

   getting_started/installation
   getting_started/running_pipeline
   getting_started/configuration

|cgr| is designed to be an easy to deploy to multiple cluster systems.
It is built on top of snakemake_ but packaged into as an easy to install python utility.
This lets us take advantage of snakemake_'s amazing workflow management system, while adding a custom library of helper tools.

.. toctree::
   :caption: QC Sub Workflows
   :name: sub-workflows
   :maxdepth: 1

   sub_workflows/entry_points
   sub_workflows/intensity_check
   sub_workflows/contamination
   sub_workflows/sample_qc
   sub_workflows/subject_qc
   sub_workflows/delivery

|cgr| is broken into 5 sub-workflows.
Each sub-workflow can be run independently, as long as all previous sub-workflows are complete.
For example, the *Sample QC* sub-workflow requires the *Entry Points* sub-workflow (and optionally the *Contamination*) to be complete.
Here we describe the various steps that each sub-workflow runs, the config options, and a summary of the generated outputs.

.. toctree::
   :caption: User Reference
   :name: module-reference
   :maxdepth: 1

   reference/file_types

.. toctree::
   :caption: Developer Reference
   :name: dev-docs
   :titlesonly:

   dev_docs/setup
   dev_docs/documentation
   dev_docs/source_code
   dev_docs/testing

.. toctree::
   :caption: API Reference
   :name: api
   :maxdepth: 1

   api/cgr_gwas_qc.rst

To Do
-----

.. todolist::
