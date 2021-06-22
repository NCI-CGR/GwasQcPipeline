.. _delivery:

Delivery Sub-workflow
=====================

**Workflow File**:
   https://github.com/NCI-CGR/GwasQcPipeline/blob/default/src/cgr_gwas_qc/workflow/sub_workflows/delivery.smk

**Config Options**: see :ref:`config-yaml` for more details

   - ``user_files.output_pattern``
   - ``workflow_params.lims_upload``

**Major Outputs**:
   - ``delivery/samples.bed``
   - ``delivery/samples.bim``
   - ``delivery/samples.fam``
   - ``delivery/SampleUsedforSubject.csv``
   - ``delivery/subjects.bed``
   - ``delivery/subjects.bim``
   - ``delivery/subjects.fam``
   - ``delivery/HWP.zip``
   - ``delivery/*AnalysisManifest*.csv``
   - ``delivery/*QC_Report*.docx``
   - ``delivery/*QC_Report*.xlsx``
   - ``files_for_lab/*all_sample_qc*.csv``
   - ``files_for_lab/*LimsUpload*.csv``
   - ``files_for_lab/*Identifiler*.csv``
   - ``files_for_lab/*KnownReplicates*.csv``
   - ``files_for_lab/*UnknownReplicates*.csv``

The delivery workflow mostly just copies files to the correct location in either ``delivery/`` or ``files_for_lab/``.
It also generates the QC report (``*.docx``) and QC report table (``*.xlsx``).
If this ``lims_upload`` options is set then it will copy ``*LimsUpload*.csv`` to the root of the project folder.
