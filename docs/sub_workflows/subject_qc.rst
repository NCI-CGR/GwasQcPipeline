.. _subject-qc:

Subject QC Sub-workflow
=======================

**Workflow File**:
   https://github.com/NCI-CGR/GwasQcPipeline/blob/default/src/cgr_gwas_qc/workflow/sub_workflows/subject_qc.smk

**Config Options**: see :ref:`config-yaml` for more details

   - ``workflow_params.minimum_pop_subject``
   - ``workflow_params.control_hwp_threshold``
   - ``software_params.ibd_pi_hat_min``
   - ``software_params.ibd_pi_hat_max``
   - ``software_params.dup_concordance_cutoff``
   - ``software_params.pi_hat_threshold``
   - ``software_params.maf_for_ibd``
   - ``software_params.maf_for_hwe``
   - ``software_params.ld_prune_r2``
   - ``software_params.autosomal_het_threshold``

**Major Outputs**:

   - ``subject_level/subject_qc.csv``
   - ``subject_level/samples.bed``
   - ``subject_level/samples.bim``
   - ``subject_level/samples.fam``
   - ``subject_level/subjects.bed``
   - ``subject_level/subjects.bim``
   - ``subject_level/subjects.fam``
   - ``subject_level/concordance.csv``
   - ``subject_level/population_qc.csv``

.. figure:: ../static/subject_qc.*
   :name: fig-subject-qc-workflow

   The subject-qc sub-workflow.
   This runs additional QC checks at the subject and population level.

Selecting Subject Representative
--------------------------------

**Major Outputs**:

   - ``subject_level/subject_qc.csv``
   - ``subject_level/samples.bed`` contains samples that are subject representatives
   - ``subject_level/samples.bim`` contains samples that are subject representatives
   - ``subject_level/samples.fam`` contains samples that are subject representatives
   - ``subject_level/subjects.bed`` renames samples to subject IDs
   - ``subject_level/subjects.bim`` renames samples to subject IDs
   - ``subject_level/subjects.fam`` renames samples to subject IDs

During :ref:`sample-qc` we build the ``sample_level/sample_qc.csv``.
In this table we have a flag ``is_subject_representative`` to indicate if a sample was selected to represent the subject.
The ``subject_level?subject_qc.csv`` is simply the ``sample_level/sample_qc.csv`` but only with subject representatives.
We then pull out these samples from ``sample_level/call_rate_2/samples.{bed,bim,fam}`` and rename them to subject IDs.

Population Level Analysis
-------------------------

**Config Options**: see :ref:`config-yaml` for more details

   - ``workflow_params.minimum_pop_subject``

**Major Outputs**:

   - ``subject_level/<population>/subjects.bed``
   - ``subject_level/<population>/subjects.bim``
   - ``subject_level/<population>/subjects.fam``

Here we split subjects into ancestral populations based on GRAF_ calls made during the :ref:`sample-qc`.

Autosomal Heterozygosity
^^^^^^^^^^^^^^^^^^^^^^^^

**Config Options**: see :ref:`config-yaml` for more details

   - ``software_params.autosomal_het_threshold``

**Major Outputs**:

   - ``subject_level/<population>/subjects.het``
   - ``subject_level/autosomal_heterozygosity_plots/<population>.png``


Here we calculate Autosomal Heterozygosity separately for each population and generate the plots used in the QC report.

Subject Relatedness
^^^^^^^^^^^^^^^^^^^

**Config Options**: see :ref:`config-yaml` for more details

   - ``software_params.ibd_pi_hat_min``
   - ``software_params.ibd_pi_hat_max``
   - ``software_params.dup_concordance_cutoff``
   - ``software_params.pi_hat_threshold``
   - ``software_params.maf_for_ibd``
   - ``software_params.ld_prune_r2``

**Major Outputs**:

   - ``subject_level/<population>/relatives.csv``
   - ``subject_level/<population>/related_subjects_to_remove.csv``

Here we estimated related (IBD) again, but separately for each population.
We then prune subjects so that no two subjects have a PI_HAT >= pi_hat_threshold.

Population Structure (PCA)
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Config Options**: see :ref:`config-yaml` for more details

   - ``software_params.maf_for_ibd``
   - ``software_params.ld_prune_r2``

**Major Outputs**:

   - ``subject_level/<population>/subjects_maf{maf_for_ibd}_ld{ld_prune_r2}.eigenvec``
   - ``subject_level/pca_plots/<population>.png``

Here we use EIGENSOFT_ fast pca to identify population structure.
We then plot a panel of pair-plots for the first size principal components.

Hardy Weinberg
^^^^^^^^^^^^^^

**Config Options**: see :ref:`config-yaml` for more details

   - ``workflow_params.control_hwp_threshold``
   - ``software_params.maf_for_hwe``

**Major Outputs**:

   - ``subject_level/<population>/controls_unrelated_maf{maf_for_hwe}_snps_autosomes.hwe``

Here we pull out only control subjects.
We then calculate Hardy Weinberg Equilibrium.
We use only controls, but cases may have SNPs that are out of HWE.

Subject QC Summary Tables
-------------------------

**Major Outputs**:

   - ``subject_level/population_qc.csv``

Finally, all results are aggregated into the poulation level QC table.

.. automodule:: cgr_gwas_qc.workflow.scripts.agg_population_qc_tables
