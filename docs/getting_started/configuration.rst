Configuration Files
===================

.. _sample-sheet:

Sample Sheet File
-----------------

A sample sheet is a CSV formatted file where each row represents a sample and columns contain various metadata.
This file must contain a column named ``Sample_ID`` which has a unique ID for each row.
You also need a column describing (1) the subject ID representing samples from an individual, (2) the expected sex of individual, and (3) the case control status of the individual.
See the `Workflow Parameters`_ section below for details about what values are allowed for subject ID, expected sex, and case control.

.. note::

   Example of the basic required structure of a sample sheet.
   You can have any other metadata as additional columns.
   The only column name that is required is ``Sample_ID``, the other columns can be named anything as long as they are referenced correctly in `Workflow Parameters`_.

   .. code-block:: CSV

      Sample_ID,Subject_ID,Expected_Sex,Case_Control
      Samp0001,Sub0001,M,Case
      Samp0002,Sub0001,M,Case
      Samp0003,Sub0002,F,Control
      Samp0004,Sub0003,F,QC

.. _lims:

CGR LIMs Manifest File
----------------------

The CGR LIMs manifest file is our internal way to distribute sample information.
If you use a `Sample Sheet File`_ then you can ignore this section.
The manifest file is an INI-like file with three section (header, manifests, data).
During ``cgr pre-flight`` we will pull out the data section.
This file should already have all of the required columns.
See the `Workflow Parameters`_ section below for details about what values are allowed for subject ID, expected sex, and case control.

.. _config-yaml:

The ``config.yml``
------------------

In :ref:`Running the Pipeline`, we created the ``config.yml`` using :ref:`cgr-config`.
This file is central to running the |cgr|.
Here we will describe the main configuration options.
The config file is broken down into different sections called "namespaces".
I will walk through each section below and end this page with a full example.

.. note::
   Required options without defaults are in bold.

.. pydantic:: cgr_gwas_qc.models.config.Config

The ``config.yml`` file contains all workflow configurations.
The top level section is automatically populated by ``cgr config`` and ``cgr pre-flight``.
You should not need to edit this section unless you want to
(1) update the version of the workflow you are using, (2) change the project name, or (3) switch human references.

.. pydantic:: cgr_gwas_qc.models.config.reference_files.ReferenceFiles

If you are on ``CGEMS/CCAD`` then the paths of the reference files are correctly populated by ``cgr config --cgems`` or ``cgr config --cgems-dev``.
If you are on another system then you most provide the correct paths and versions of these files.
The ``illumina_manifest_file`` is provided by Illumina.
The ``thousand_genome_vcf`` and ``thousand_genome_tbi`` files can be downloaded from the 1000 Genome's website.
The ``illumina_cluster_file`` is the EGT file used to generate GTCs; this file is not required and is only referenced in the QC report if present.

1000 Genomes reference files download links
   * **hg37**: `VCF <http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz>`__,
     `TBI <http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz.tbi>`__

   * **hg38**: `VCF <http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz>`__,
     `TBI <http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz.tbi>`__

.. warning::

   The genome build needs to match ``illumina_manifest_file``, ``thousand_genome_*``, and ``genome_build`` above.
   Currently we only support ``hg37`` and ``hg38``.
   Though in reality as long as ``illumina_manifest_file`` and ``thousand_genome_*`` are from the same build then everything should work.

.. pydantic:: cgr_gwas_qc.models.config.user_files.UserFiles

The ``user_files`` section has a number of different mutually exclusive configurations depending on the starting file types.
By default we assume you will be starting with IDAT and GTC files.
Though we also accept aggregated PED/MAP and aggregated BED/BIM/FAM files.
In the example, ``{Project}`` and ``{Sample_ID}`` will be filled by values from ``Project`` and ``Sample_ID`` columns in ``cgr_sample_sheet.csv``.

If the ``gtc_pattern`` is given, then this will trigger the GTC entry point.
We will convert each sample's GTC to a PED/MAP and then aggregate all samples and convert to a single BED/BIM/FAM at ``sample_level/samples.{bed,bim,fam}``.

If the PED/MAP files are given, then this will trigger the PED/MAP entry point.
Which will convert these files to a single BED/BIM/FAM at ``sample_level/samples.{bed,bim,fam}``.

If the BED/BIM/FAM files are given, then this will trigger the BED/BIM/FAM entry point.
Which will create a symbolic link from your BED/BIM/FAM to ``sample_level/samples.{bed,bim,fam}``.

.. pydantic:: cgr_gwas_qc.models.config.software_params.SoftwareParams

.. pydantic:: cgr_gwas_qc.models.config.workflow_params.WorkflowParams


Sample IDs to Remove
^^^^^^^^^^^^^^^^^^^^

This is an optional section where you can list ``Sample_ID`` that you do not want to include in the QC run.
These samples will be indicated as ``is_user_exclusion = True`` in the sample level QC table.

.. code-block:: yaml

   Sample_IDs_to_remove:
      - Sample0001

Full Example
^^^^^^^^^^^^

.. code-block:: yaml

   pipeline_version: v1.0.0
   project_name: SR0001-001_1_0000000
   sample_sheet: /path/to/manifest/file/SR0001-001_1_AnalysisManifest_0000000.csv
   genome_build: hg37
   snp_array: GSAMD-24v1-0
   num_samples: 6336
   num_snps: 700078
   reference_files:
      illumina_manifest_file: /path/to/bpm/file/GSAMD-24v1-0_20011747_A1.bpm
      thousand_genome_vcf: /path/to/thousand/genome/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz
      thousand_genome_tbi: /path/to/thousand/genome/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz.tbi
   user_files:
      output_pattern: '{prefix}/{file_type}.{ext}'
      idat_pattern:
         red: /expample/pattern/wildcards/are/columns/in/sample_sheet_file/{Project}/{Sample_ID}_Red.idat
         green: /expample/pattern/wildcards/are/columns/in/sample_sheet_file/{Project}/{Sample_ID}_Grn.idat
      gtc_pattern: /expample/pattern/wildcards/are/columns/in/sample_sheet_file/{Project}/{Sample_ID}.gtc
   software_params:
      sample_call_rate_1: 0.8
      snp_call_rate_1: 0.8
      sample_call_rate_2: 0.95
      snp_call_rate_2: 0.95
      ld_prune_r2: 0.1
      maf_for_ibd: 0.2
      maf_for_hwe: 0.05
      ibd_pi_hat_min: 0.05
      ibd_pi_hat_max: 1.0
      dup_concordance_cutoff: 0.95
      intensity_threshold: 6000
      contam_threshold: 0.1
      contam_population: AF
      pi_hat_threshold: 0.2
      autosomal_het_threshold: 0.1
      strand: top
   workflow_params:
      subject_id_column: Group_By
      expected_sex_column: Expected_Sex
      case_control_column: Case/Control_Status
      remove_contam: true
      remove_rep_discordant: true
      minimum_pop_subjects: 50
      control_hwp_threshold: 50
      lims_upload: false
   Sample_IDs_to_remove:
      - Sample0001
