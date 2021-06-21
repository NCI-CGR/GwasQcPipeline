.. _sample-qc:

Sample QC Sub-workflow
======================

**Workflow File**:
   https://github.com/NCI-CGR/GwasQcPipeline/blob/default/src/cgr_gwas_qc/workflow/sub_workflows/sample_qc.smk

**Config Options**: see :ref:`config-yaml` for more details

   - ``reference_files.thousand_genome_vcf``
   - ``reference_files.thousand_genome_tbi``
   - ``software_params.sample_call_rate_1``
   - ``software_params.snp_call_rate_1``
   - ``software_params.sample_call_rate_2``
   - ``software_params.snp_call_rate_2``
   - ``software_params.intensity_threshold``
   - ``software_params.contam_threshold``
   - ``software_params.ibd_pi_hat_min``
   - ``software_params.ibd_pi_hat_max``
   - ``software_params.maf_for_ibd``
   - ``software_params.ld_prune_r2``
   - ``software_params.dup_concordance_cutoff``
   - ``software_params.pi_hat_threshold``
   - ``workflow_params.remove_contam``
   - ``workflow_params.remove_rep_discordant``

**Major Outputs**:

   - ``sample_level/sample_qc.csv``
   - ``sample_level/snp_qc.csv``
   - ``sample_level/concordance/KnownReplicates.csv``
   - ``sample_level/concordance/InternalQcKnown.csv``
   - ``sample_level/concordance/StudySampleKnown.csv``
   - ``sample_level/concordance/UnknownReplicates.csv``
   - ``sample_level/summary_stats.txt``
   - ``sample_level/call_rate.png``
   - ``sample_level/chrx_inbreeding.png``
   - ``sample_level/ancestry.png``

.. figure:: ../static/sample_qc.*
   :name: fig-sample-qc-workflow

   The sample-qc sub-workflow.
   This runs all major QC checks at the sample level.


Call Rate Filters
-----------------

**Config Options**: see :ref:`config-yaml` for more details

   - ``software_params.sample_call_rate_1``
   - ``software_params.snp_call_rate_1``
   - ``software_params.sample_call_rate_2``
   - ``software_params.snp_call_rate_2``

**Major Outputs**:

   - ``sample_level/call_rate_1/samples.{bed,bim,fam}`` data set after level 1 sample/snp call rate filters.
   - ``sample_level/call_rate_2/samples.{bed,bim,fam}`` data set after level 2 sample/snp call rate filters.

First we run a series of call rate (a.k.a. completion) filters.
Here we remove samples then SNPs which have high missingness in across samples.


Contamination Summary
---------------------

**Config Options**: see :ref:`config-yaml` for more details

   - ``software_params.intensity_threshold``
   - ``software_params.contam_threshold``

**Major Outputs**:

   - ``sample_level/contamination/summary.csv`` aggregate summary table of contamination results.

If IDAT/GTC files are present and we ran :ref:`contamination` then we will create a summary table.
For samples that have a low call rate and a median idat intensity < the intensity_threshold we will set %Mix to NaN.
We then will flag a sample ``is_contaminated`` if the %Mix is >= contam_threshold.

Replicate Concordance
---------------------

**Config Options**: see :ref:`config-yaml` for more details

   - ``software_params.ibd_pi_hat_min``
   - ``software_params.ibd_pi_hat_max``
   - ``software_params.maf_for_ibd``
   - ``software_params.ld_prune_r2``
   - ``software_params.dup_concordance_cutoff``
   - ``software_params.pi_hat_threshold``
   - ``workflow_params.remove_contam``
   - ``workflow_params.remove_rep_discordant``

**Major Outputs**:

   - ``sample_level/concordance/summary.csv``

It is common to include multiple replicates for a subset of samples as a QC metric.
Here we check these replicates and make sure that they are highly concordant.
We report 3 methods for checking relatedness/concordance (PLINK IBD, GRAF, and KING).
We mostly rely on IBD results but report the others in our summary tables.
For all three methods, we do all pairwise comparisons to determine which samples are related/replicated.
We flag two samples as concordant (``is_ge_concordance``) if concordance (proportion IBD2) > dup_concordance_cutoff.
We flag two samples as related (``is_ge_pi_hat``) if PLINK's PI_HAT is >= pi_hat_threshold.

Ancestry
--------

**Major Outputs**:

   - ``sample_level/ancestry/graf_ancestry.txt``

Finally, we estimate ancestry using GRAF_.
Here we categorize each sample into 1 of 9 possible groups (see GRAF_ for more details).
We also report the percent mixture of the 3 major continental population (African, East Asian, European).

Sample QC Summary
-----------------

**Major Outputs**:

   - ``sample_level/sample_qc.csv``
   - ``sample_level/snp_qc.csv``
   - ``sample_level/summary_stats.txt``
   - ``sample_level/call_rate.png``
   - ``sample_level/chrx_inbreeding.png``
   - ``sample_level/ancestry.png``

Finally, we aggregate all of the QC results.
We generate two summary tables described below:

.. _sample-qc-table:

.. automodule:: cgr_gwas_qc.workflow.scripts.sample_qc_table

.. _snp-qc-table:

.. automodule:: cgr_gwas_qc.workflow.scripts.snp_qc_table
