"""This module contains sample level filters."""

if (
    cfg.config.user_files.idat_pattern
    and cfg.config.user_files.gtc_pattern
    and cfg.config.workflow_params.remove_contam
):

    ################################################################################
    # Contaminated Samples
    ################################################################################

    rule idat_intensity:
        """Calculate the median intensity overall intensity (Red + Green).

        .. warning::
            This is a submission hot-spot creating one job per sample. Each output file contains a
            single number, the median intensity.
        """
        input:
            red=cfg.config.user_files.idat_pattern.red,
            green=cfg.config.user_files.idat_pattern.green,
        output:
            temp("idat_intensity/{Sample_ID}.{SentrixBarcode_A}.{SentrixPosition_A}.txt"),
        envmodules:
            cfg.envmodules("r"),
        conda:
            "../conda/illuminaio.yml"
        script:
            "../scripts/median_idat_intensity.R"

    rule combine_idat_intensity:
        """Aggregates sample level median intensity values into a single table."""
        input:
            cfg.expand(rules.idat_intensity.output[0]),
        output:
            "all_sample_idat_intensity/idat_intensity.csv",
        run:
            with open(output[0], "w") as out:
                out.write("SampId,ChipId,MedianIntensity\n")
                for i in input:
                    pth = Path(i)
                    sample_id, barcode, position = pth.stem.split(".")
                    median_intensity = pth.read_text().strip()
                    out.write(f"{sample_id},{barcode}_{position},{median_intensity}\n")

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

    rule combine_sample_level_contam_test:
        """Aggregate sample contamination scores.

        Aggregates sample level contamination scores into a single file (each row is a sample). The
        script sets `%Mix` to `NA` if the intensity is below the threshold and the file is not in the
        `imiss3` file.
        """
        input:
            contam=cfg.expand("one_samp_b_1000g_contam/{Sample_ID}.contam.out"),
            intens="all_sample_idat_intensity/idat_intensity.csv",
            imiss3="plink_filter_call_rate_2/samples_filter2.imiss",
        params:
            intensThresh=6000,
        output:
            "all_contam/contam.csv",
        run:
            crDict3 = makeCallRateDict(input.imiss3)
            intensDict = {}
            with open(input.intens) as f:
                head = f.readline()
                line = f.readline()
                while line != "":
                    (samp, chipId, intensity) = line.rstrip().split(",")
                    intensDict[samp] = float(intensity)
                    line = f.readline()
            with open(output[0], "w") as out:
                out.write("ID,%Mix,LLK,LLK0\n")
                for i in input.contam:
                    samp = os.path.basename(i).split(".")[0]
                    intens = intensDict[samp]
                    with open(i) as f:
                        head = f.readline()
                        while "%Mix" not in head and head != "":
                            head = f.readline()
                        if head == "":
                            print("strange file format: " + i)
                            sys.exit(1)
                        head = f.readline()
                        line = f.readline()
                        line_list = line.split()
                        line_list[0] = samp
                        if intens < params.intensThresh and not crDict3.get(samp):
                            line_list[1] = "NA"
                        out.write(",".join(line_list) + "\n")
