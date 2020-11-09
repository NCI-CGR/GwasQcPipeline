"""This module contains sample level filters."""

if cfg.config.user_files.gtc_patterns and cfg.config.workflow_params.remove_contam:

    ################################################################################
    # Contaminated Samples
    ################################################################################
    rule gtc_to_adpc:
        """Converts a sample's GTC/BPM to an Illumina ADPC.BIN.

        This is the format required by ``verifyIDintensity``. The script also runs some sanity checks
        (intensities and normalized intensities > 0; genotypes one of {0, 1, 2, 3}) while processing
        each file.

        .. warning::
            This is a submission hot spot creating 1 job per sample.
        """
        input:
            gtc=lambda wc: cfg.expand(
                cfg.config.user_files.gtc_pattern, query=f"Sample_ID == '{wc.Sample_ID}'",
            )[0],
            bpm=cfg.config.reference_files.illumina_manifest_file,
        output:
            adpc=temp("contam/{Sample_ID}.adpc.bin"),
            snp_count="contam/{Sample_ID}.adpc.bin.numSnps.txt",
        script:
            "../scripts/gtc2adpc.py"

    rule pull_1KG_allele_b_freq:
        """Pulls the population level allele frequencies from the 1KG project.

        ``verifyIDintensity`` requires population level allele frequencies for its model. Here we use
        a custom script to pull out the allele B frequencies (ABF) from the 1000 genomes project
        (1KG). To do this we take each marker from the manifest file (BPM) and pull out ABF in the
        1KG ``.vcf`` from the ``INFO`` column. The script allows pulling out allele frequencies for
        different super populations but defaults to ``AF`` which ignores super population.
        """
        input:
            bpm=cfg.config.reference_files.illumina_manifest_file,
        params:
            population=cfg.config.software_params.contam_population,
        output:
            (
                cfg.config.reference_files.illumina_manifest_file.stem
                + cfg.config.software_params.contam_population
                + "abf.txt",
            ),
        script:
            "bpm2abf.py"

    rule sample_level_contam_test:
        """Find contaminated samples using allele intensities.

        Uses ``verifyIDintensity`` to find samples with allele intensities that deviate from the
        population.

        .. warning::
            This is a submission hot spot creating 1 job per sample.

        .. note::
            Here we are running ``verifyIDintensity`` in single sample mode. This software also has a
            multi-sample mode which may be faster and give better estimates. The problem with
            multi-sample mode is that it only works when you have a "large" number of samples.
        """
        input:
            adpc=rules.gtc_to_adpc.output.adpc,
            abf=rules.pull_1KG_allele_b_freq.output[0],
        params:
            snps=numSNPs,
        output:
            temp("one_samp_b_1000g_contam/{Sample_ID}.contam.out"),
        conda:
            "../conda/verifyidintensity.yml"
        shell:
            "verifyIDintensity -m {params.snps} -n 1 -b {input.abf} -v -p -i {input.adpc} > {output}"
