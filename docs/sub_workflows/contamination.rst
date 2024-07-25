.. _contamination:

Contamination Sub-workflow
==========================

**Workflow File**:
   https://github.com/NCI-CGR/GwasQcPipeline/blob/default/src/cgr_gwas_qc/workflow/sub_workflows/contamination.smk

**Config Options**: see :ref:`config-yaml` for more details

   - ``reference_files.thousand_genome_vcf``
   - ``reference_files.thousand_genome_tbi``
   - ``user_files.gtc_pattern``
   - ``user_files.idat_pattern``
   - ``software_params.contam_population``

**Major Outputs**:

   - ``sample_level/<BPM Prefix>.<software_params.contam_population>.abf.txt`` B allele frequencies from the 1000 genomes.
   - ``sample_level/contamination/median_idat_intensity.csv`` aggregated table of median IDAT intensities.
   - ``sample_level/contamination/verifyIDintensity.csv`` aggregated table of contamination scores.


.. figure:: ../static/contamination.png
   :name: fig-contamination-workflow

   The contamination sub-workflow.
   This workflow will estimate contamination using verifyIDintensity on each sample individually.
   It requires that you have GTC/IDAT files.
   It first pulls B-allele frequencies from the 1000 Genomes VCF file.
   It then estimate contamination for each sample and aggregates these results.
   Finally, it also estimates the per sample median IDAT intensity, which is used to filter contamination results in the :ref:`sample-qc`
