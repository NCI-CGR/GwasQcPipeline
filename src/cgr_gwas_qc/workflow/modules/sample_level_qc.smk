import pandas as pd


include: cfg.modules("common.smk")


################################################################################
# Sample Contamination
#
# Requires:
#   - IDAT Files (red and green)
#   - GTC Files
################################################################################
if (
    cfg.config.user_files.idat_pattern
    and cfg.config.user_files.gtc_pattern
    and cfg.config.workflow_params.remove_contam
):

    rule per_sample_median_idat_intensity:
        """Calculate median intensity of Red and Green channels."""
        input:
            red=lambda wc: cfg.expand(
                cfg.config.user_files.idat_pattern.red, query=f"Sample_ID == '{wc.Sample_ID}'",
            ),
            green=lambda wc: cfg.expand(
                cfg.config.user_files.idat_pattern.green, query=f"Sample_ID == '{wc.Sample_ID}'",
            ),
        output:
            temp(
                "sample_level/per_sample_median_idat_intensity/{Sample_ID}.{SentrixBarcode_A}.{SentrixPosition_A}.csv"
            ),
        resources:
            mem_gb=1,
        group:
            "per_sample_median_idat_intensity"
        envmodules:
            cfg.envmodules("r"),
        conda:
            cfg.conda("illuminaio.yml")
        script:
            "../scripts/median_idat_intensity.R"

    rule agg_median_idat_intensity:
        input:
            cfg.expand(rules.per_sample_median_idat_intensity.output[0]),
        output:
            "sample_level/median_idat_intensity.csv",
        run:
            pd.concat([pd.read_csv(file_name) for file_name in input]).to_csv(
                output[0], index=False
            )

    rule per_sample_gtc_to_adpc:
        """Converts a sample's GTC/BPM to an Illumina ADPC.BIN.

        This is the format required by ``verifyIDintensity``. The script also
        runs some sanity checks (intensities and normalized intensities > 0;
        genotypes are one of {0, 1, 2, 3}) while processing each file.

        .. warning::
            This is a submission hot spot creating 1 job per sample.
        """
        input:
            gtc_file=lambda wc: cfg.expand(
                cfg.config.user_files.gtc_pattern, query=f"Sample_ID == '{wc.Sample_ID}'",
            )[0],
            bpm_file=cfg.config.reference_files.illumina_manifest_file,
        output:
            temp("sample_level/per_sample_adpc/{Sample_ID}.adpc.bin"),
        resources:
            mem_gb=1,
        group:
            "per_sample_gtc_to_adpc"
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
            bpm_file=cfg.config.reference_files.illumina_manifest_file,
            vcf_file=cfg.config.reference_files.thousand_genome_vcf,
        params:
            population=cfg.config.software_params.contam_population,
        output:
            abf_file="sample_level/{}.{}.abf.txt".format(
                cfg.config.reference_files.illumina_manifest_file.stem,
                cfg.config.software_params.contam_population,
            ),
        script:
            "../scripts/bpm2abf.py"

    rule per_sample_verifyIDintensity_contamination:
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
            adpc=rules.per_sample_gtc_to_adpc.output[0],
            abf=rules.pull_1KG_allele_b_freq.output.abf_file,
        params:
            snps=cfg.config.num_snps,
        output:
            temp("sample_level/per_sample_contamination_test/{Sample_ID}.contam.out"),
        resources:
            mem_gb=1,
        group:
            "per_sample_verifyIDintensity_contamination"
        conda:
            cfg.conda("verifyidintensity.yml")
        shell:
            "verifyIDintensity -m {params.snps} -n 1 -b {input.abf} -v -p -i {input.adpc} > {output}"

    rule agg_verifyIDintensity_contamination:
        """Aggregate sample contamination scores.

        Aggregates sample level contamination scores into a single file (each
        row is a sample). The script sets ``%Mix`` to ``NA`` if the intensity
        is below the threshold and the file is not in the ``imiss`` file.
        """
        input:
            contamination=cfg.expand(rules.per_sample_verifyIDintensity_contamination.output),
            median_idat_intensity=rules.agg_median_idat_intensity.output[0],
            imiss="sample_level/call_rate_2/samples.imiss",
        params:
            intensity_threshold=cfg.config.software_params.intensity_threshold,
        output:
            "sample_level/contamination/verifyIDintensity_contamination.csv",
        script:
            "../scripts/agg_contamination_test.py"


################################################################################
# GRAF Conversions
################################################################################
rule update_snps_to_1kg_rsID:
    """Update SNP IDs to rsID from the 1KG project.

    Update study marker IDs to correspond with GRAF's fingerprints which are
    based on the 1KG project rsIDs.
    """
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
        vcf=cfg.config.reference_files.thousand_genome_vcf,
    output:
        bed="{prefix}_1kg_rsID.bed",
        bim="{prefix}_1kg_rsID.bim",
        fam="{prefix}_1kg_rsID.fam",
        id_map="{prefix}_1kg_rsID.csv",
    script:
        "../scripts/update_snps_to_1kg_rsID.py"


rule graf_extract_fingerprint_snps:
    """Extract fingerprint SNPs from a PLINK data and convert to a GRAF."""
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
    output:
        "{prefix}.fpg",
    log:
        "{prefix}.fpg.log",
    conda:
        cfg.conda("graf.yml")
    shell:
        "graf "
        "-exfp {wildcards.prefix} "
        "-out {output[0]} "
        "-type 4 "
        "> {log} 2>&1 "
        "|| exit_code=$?; if [ $exit_code -ne 1 ]; then exit $exit_code; fi" # GRAF returns an exit code of 1, this captures it so snakemake will actually run.


################################################################################
# Sample/Replicate Concordance
################################################################################
rule sample_concordance_plink:
    """Summarize sample concordance.

    Splits the concordance table into known and unknown replicates.
    """
    input:
        "sample_level/call_rate_2/samples_maf{maf}_ld{ld}_pruned.concordance.csv".format(
            maf=cfg.config.software_params.maf_for_ibd, ld=cfg.config.software_params.ld_prune_r2,
        ),
    output:
        "sample_level/concordance/plink.csv",
    shell:
        "cp {input[0]} {output[0]}"


rule sample_concordance_graf:
    """Estimate relatedness among samples.

    Outputs a table with pairwise samples and their genotypic relationship.
    GRAF uses a set of 10K pre-screened SNPs, so we can directly use Call
    Rate 2 filtered samples.
    """
    input:
        fpg="sample_level/call_rate_2/samples_1kg_rsID.fpg",
    output:
        "sample_level/concordance/graf.tsv",
    conda:
        cfg.conda("graf.yml")
    log:
        "sample_level/concordance/graf.log",
    shell:
        "graf "
        "-geno {input.fpg} "
        "-type 4 "
        "-out {output[0]} "
        "> {log} 2>&1 "
        "|| exit_code=$?; if [ $exit_code -ne 1 ]; then exit $exit_code; fi" # GRAF returns an exit code of 1, this captures it so snakemake will actually run.


rule sample_concordance_king:
    input:
        bed="sample_level/call_rate_2/samples.bed",
        bim="sample_level/call_rate_2/samples.bim",
        fam="sample_level/call_rate_2/samples.fam",
    output:
        within_family="sample_level/concordance/king.kin",
        between_family="sample_level/concordance/king.kin0",
        within_family_X="sample_level/concordance/kingX.kin",
        between_family_X="sample_level/concordance/kingX.kin0",
    params:
        out_prefix="sample_level/concordance/king",
    conda:
        cfg.conda("king.yml")
    log:
        "sample_level/concordance/king.log",
    shell:
        "king -b {input.bed} --kinship --prefix {params.out_prefix} > {log} 2>&1 "
        "&& touch {output.within_family} "
        "&& touch {output.between_family} "
        "&& touch {output.within_family_X} "
        "&& touch {output.between_family_X} " # king does not always output all files
         # adding these touch statements so it plays nice
         # with snakemake


rule sample_concordance:
    input:
        sample_sheet_csv="cgr_sample_sheet.csv",
        plink_file=rules.sample_concordance_plink.output[0],
        graf_file=rules.sample_concordance_graf.output[0],
        king_file=rules.sample_concordance_king.output.between_family,
    output:
        "sample_level/concordance/summary.csv",
    script:
        "../scripts/sample_concordance.py"


rule split_sample_concordance:
    input:
        rules.sample_concordance.output[0],
    output:
        known_csv="sample_level/concordance/KnownReplicates.csv",
        known_qc_csv="sample_level/concordance/InternalQcKnown.csv",
        known_study_csv="sample_level/concordance/StudySampleKnown.csv",
        unknown_csv="sample_level/concordance/UnknownReplicates.csv",
    script:
        "../scripts/split_sample_concordance.py"


################################################################################
# Sample Level Ancestry
################################################################################
rule graf_ancestry:
    """Estimate ancestry for each sample."""
    input:
        fpg="sample_level/call_rate_2/samples_1kg_rsID.fpg",
    output:
        "sample_level/ancestry/graf_populations.txt",
    log:
        "sample_level/ancestry/graf_populations.log",
    envmodules:
        cfg.envmodules("graf"),
    conda:
        cfg.conda("graf.yml")
    shell:
        "graf "
        "-geno {input.fpg} "
        "-pop {output[0]} "
        "> {log} 2>&1 "
        "|| exit_code=$?; if [ $exit_code -ne 1 ]; then exit $exit_code; fi" # GRAF returns an exit code of 1, this captures it so snakemake will actually run.


rule graf_ancestry_table:
    """Create summary table with ancestry calls."""
    input:
        rules.graf_ancestry.output[0],
    output:
        "sample_level/ancestry/graf_ancestry.txt",
    envmodules:
        cfg.envmodules("graf"),
    conda:
        cfg.conda("graf.yml")
    shell:
        "PlotPopulations.pl {input[0]} {output[0]} "
