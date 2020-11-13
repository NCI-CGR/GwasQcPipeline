"""This module contains sample level filters."""

if (
    cfg.config.user_files.idat_pattern
    and cfg.config.user_files.gtc_pattern
    and cfg.config.workflow_params.remove_contam
):
    # Remove Samples with high contamination.
    #
    # If we have IDAT and GTC files we can estimated sample contamination by
    # modeling SNP intensities. Here we remove samples that have a high
    # estimated contamination rate (config.software_params.contam_threshold)
    # [default: 0.2].
    #
    # For contamination estimation we ignore samples whose median intensity is
    # <6,000. We also use B-allele frequencies from the 1,000 genomes project
    # as part of the modeling (config.software_params.contam_population)
    # [default: AF].
    #
    # .. warning::
    #     The default contamination rate of 20% seems really high.
    #

    rule median_idat_intensity:
        """Calculate median intensity overall intensity of Red + Green channels.

        .. warning::
            This is a submission hot-spot creating one job per sample. Each output file contains a
            single number, the median intensity.
        """
        input:
            red=cfg.config.user_files.idat_pattern.red,
            green=cfg.config.user_files.idat_pattern.green,
        output:
            temp(
                "sample_filters/median_idat_intensity/{Sample_ID}.{SentrixBarcode_A}.{SentrixPosition_A}.txt"
            ),
        envmodules:
            cfg.envmodules("r"),
        conda:
            "../conda/illuminaio.yml"
        script:
            "../scripts/median_idat_intensity.R"

    rule agg_median_idat_intensity:
        """Aggregates sample level median intensity values into a single table."""
        input:
            cfg.expand(rules.median_idat_intensity.output[0]),
        output:
            "sample_filters/agg_median_idat_intensity.csv",
        run:
            with open(output[0], "w") as out:
                out.write("Sample_ID,Chip_ID,median_intensity\n")
                for i in input:
                    pth = Path(i)
                    sample_id, barcode, position = pth.stem.split(".")
                    median_intensity = pth.read_text().strip()
                    out.write(f"{sample_id},{barcode}_{position},{median_intensity}\n")

    rule convert_gtc_to_illumina_adpc:
        """Converts a sample's GTC/BPM to an Illumina ADPC.BIN.

        This is the format required by ``verifyIDintensity``. The script also
        runs some sanity checks (intensities and normalized intensities > 0;
        genotypes are one of {0, 1, 2, 3}) while processing each file.

        .. warning::
            This is a submission hot spot creating 1 job per sample.
        """
        input:
            gtc=lambda wc: cfg.expand(
                cfg.config.user_files.gtc_pattern, query=f"Sample_ID == '{wc.Sample_ID}'",
            )[0],
            bpm=cfg.config.reference_files.illumina_manifest_file,
        output:
            adpc=temp("sample_filters/convert_gtc_to_illumina_adpc/{Sample_ID}.adpc.bin"),
            snp_count=(
                "sample_filters/convert_gtc_to_illumina_adpc/{Sample_ID}.adpc.bin.numSnps.txt"
            ),
        script:
            "../scripts/gtc2adpc.py"

    rule pull_1KG_allele_b_freq:
        """Pulls the population level allele frequencies from the 1KG project.

        ``verifyIDintensity`` requires population level allele frequencies
        for its model. Here we use a custom script to pull out the allele B
        frequencies (ABF) from the 1000 genomes project (1KG). To do this we
        take each marker from the manifest file (BPM) and pull out ABF in the
        1KG ``.vcf`` from the ``INFO`` column. The script allows pulling out
        allele frequencies for different super populations but defaults to
        ``AF`` which ignores super population.
        """
        input:
            bpm=cfg.config.reference_files.illumina_manifest_file,
        params:
            population=cfg.config.software_params.contam_population,
        output:
            "sample_filters/{}.{}.abf.txt".format(
                cfg.config.reference_files.illumina_manifest_file.stem,
                cfg.config.software_params.contam_population,
            ),
        script:
            "bpm2abf.py"

    rule contamination_test:
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
            adpc=rules.convert_gtc_to_illumina_adpc.output.adpc,
            abf=rules.pull_1KG_allele_b_freq.output[0],
        params:
            snps=numSNPs,
        output:
            temp("sample_filters/contamination_test/{Sample_ID}.contam.out"),
        conda:
            "../conda/verifyidintensity.yml"
        shell:
            "verifyIDintensity -m {params.snps} -n 1 -b {input.abf} -v -p -i {input.adpc} > {output}"

    rule agg_contamination_test:
        """Aggregate sample contamination scores.

        Aggregates sample level contamination scores into a single file (each row is a sample). The
        script sets `%Mix` to `NA` if the intensity is below the threshold and the file is not in the
        `imiss3` file.
        """
        input:
            contamination=cfg.expand(rules.contamination_test.output),
            median_idat_intensity=rules.agg_median_idat_intensity.output[0],
            imiss="plink_filter_call_rate_2/samples.imiss",
        params:
            intensity_threshold=cfg.config.software_params.intensity_threshold,
        output:
            "sample_filters/agg_contamination_test.csv",
        script:
            "../scripts/agg_contamination_test.py"
