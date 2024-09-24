.. _entry-points:

Entry Points Sub-workflow
=========================

**Workflow File**:
   https://github.com/NCI-CGR/GwasQcPipeline/blob/default/src/cgr_gwas_qc/workflow/sub_workflows/entry_points.smk

**Config Options**: see :ref:`config-yaml` for more details

   - ``user_files.gtc_pattern``
   - ``user_files.idat_pattern``
   - ``user_files.ped``
   - ``user_files.map``
   - ``user_files.bed``
   - ``user_files.bim``
   - ``user_files.fam``
   - ``user_files.bcf``

**Major Outputs**:

   - ``sample_level/samples.bed``
   - ``sample_level/samples.bim``
   - ``sample_level/samples.fam``

There are four paths we can take to create these files:

1. If GTC files are provided using ``user_files.gtc_pattern`` then we will

.. figure:: ../static/entry-points_gtc.png
   :name: fig-entry-points-gtc

   The GTC entry-point workflow.
   If per sample GTC files are provided, then we will convert these files to the PED/MAP format and merge them together.

2. If an aggregated PED/MAP is provided using ``user_files.ped`` and ``user_files.map`` then we will convert the PED/MAP to BED/BIM/FAM.
3. If an aggregated BED/BIM/FAM is provided using ``user_files.bed``, ``user_files.bim``, ``user_files.fam`` then we will create a symbolic link.
4. If an aggregated BCF file is provided using ``user_files.bcf`` then we will convert the BCF to BED/BIM/FAM.
