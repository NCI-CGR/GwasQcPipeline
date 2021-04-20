Command Line Interface
======================

.. warning::
    ``--ibd_pi_hat_cutoff`` is not used anywhere.

.. code-block:: bash

    ./GwasQcWrapper.py -h
    usage: GwasQcWrapper.py [-h] [-p PATH_TO_PLINK_FILE] [-d DIRECTORY_FOR_OUTPUT]
                            [--snp_cr_1 SNP_CR_1] [--samp_cr_1 SAMP_CR_1]
                            [--snp_cr_2 SNP_CR_2] [--samp_cr_2 SAMP_CR_2]
                            [--ld_prune_r2 LD_PRUNE_R2]
                            [--maf_for_ibd MAF_FOR_IBD] -s SAMPLE_SHEET
                            [--subject_id_column SUBJECT_ID_COLUMN]
                            [--ibd_pi_hat_cutoff IBD_PI_HAT_CUTOFF]
                            [--dup_concordance_cutoff DUP_CONCORDANCE_CUTOFF]
                            [--lims_output_dir LIMS_OUTPUT_DIR]
                            [--contam_threshold CONTAM_THRESHOLD]
                            [-i ILLUMINA_MANIFEST_FILE] [-a ADPC_FILE]
                            [-g GTC_DIR]
                            [--expected_sex_column EXPECTED_SEX_COLUMN]
                            [-q QUEUE] [-u]

    optional arguments:
    -h, --help            show this help message and exit
    -i ILLUMINA_MANIFEST_FILE, --illumina_manifest_file ILLUMINA_MANIFEST_FILE
                            Full path to illimina .bpm manifest file. Required for
                            gtc files.
    -a ADPC_FILE, --adpc_file ADPC_FILE
                            Full path to adpc.bin file. Required for PLINK input.
    -g GTC_DIR, --gtc_dir GTC_DIR
                            Full path to gtc directory to use instead of project
                            directory, which is the default. Will recursively find
                            gtc files in this directory.
    -u, --unlock_snakemake
                            OPTIONAL. If pipeline was killed unexpectedly you may
                            need this flag to rerun

    Required Arguments:
    -d DIRECTORY_FOR_OUTPUT, --directory_for_output DIRECTORY_FOR_OUTPUT
                            REQUIRED. Full path to the base directory for the Gwas
                            QC pipeline output. Defaults to
                            /DCEG/CGF/GWAS/Scans/GSA_Lab_QC/SR#/builds/QC_v#_date
    -s SAMPLE_SHEET, --sample_sheet SAMPLE_SHEET
                            Full path to illimina style sample sheet csv file.
    --subject_id_column SUBJECT_ID_COLUMN
                            Name of column in sample sheet that corresponds to
                            subject ID to use.
    --expected_sex_column EXPECTED_SEX_COLUMN
                            Name of column in sample sheet that corresponds to
                            expected sex of sample.

    Exactly one of these arguments is required:
    -p PATH_TO_PLINK_FILE, --path_to_plink_file PATH_TO_PLINK_FILE
                            Full path to either PLINK ped or bed to use as input.
                            need either this or gtc file project directory -g

    Required arguments with default settings:
    --snp_cr_1 SNP_CR_1   REQUIRED. SNP call rate filter 1. default= 0.80
    --samp_cr_1 SAMP_CR_1
                            REQUIRED. Sample call rate filter 1. default= 0.80
    --snp_cr_2 SNP_CR_2   REQUIRED. SNP call rate filter 2. default= 0.95
    --samp_cr_2 SAMP_CR_2
                            REQUIRED. Sample call rate filter 2. default= 0.95
    --ld_prune_r2 LD_PRUNE_R2
                            REQUIRED. r-squared cutoff for ld pruning of SNPs to
                            use for IBD and concordance. default= 0.10
    --maf_for_ibd MAF_FOR_IBD
                            REQUIRED. MAF cutoff of SNPs to use for IBD and
                            concordance. default= 0.20
    --ibd_pi_hat_cutoff IBD_PI_HAT_CUTOFF
                            REQUIRED. PI_HAT cutoff to call samples replicates.
                            default= 0.95
    --dup_concordance_cutoff DUP_CONCORDANCE_CUTOFF
                            REQUIRED. SNP concordance cutoff to call samples
                            replicates. default= 0.95
    --lims_output_dir LIMS_OUTPUT_DIR
                            Directory to copy QC file to upload to LIMS
    --contam_threshold CONTAM_THRESHOLD
                            REQUIRED. Cutoff to call a sample contaminated.
                            default= 0.10
    -q QUEUE, --queue QUEUE
                            OPTIONAL. Queue on cgemsiii to use to submit jobs.
                            Defaults to all of the seq queues and all.q if not
                            supplied. default="all.q,seq-alignment.q,seq-calling.q
                            ,seq-calling2.q,seq-gvcf.q"
