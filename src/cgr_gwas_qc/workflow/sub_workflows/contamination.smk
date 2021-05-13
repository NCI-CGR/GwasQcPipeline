import pandas as pd

from cgr_gwas_qc import load_config

cfg = load_config()


################################################################################
# Contamination Targets
################################################################################
rule all_contamination:
    input:
        "sample_level/contamination/median_idat_intensity.csv",
        "sample_level/contamination/verifyIDintensity.csv",


################################################################################
# Imports
################################################################################
module thousand_genomes:
    snakefile:
        cfg.modules("thousand_genomes")
    config:
        {}


use rule pull_b_allele_freq_from_1kg from thousand_genomes


################################################################################
# Workflow Rules
################################################################################
# -------------------------------------------------------------------------------
# Calculate Median Idat Intensity
# -------------------------------------------------------------------------------
rule per_sample_median_idat_intensity:
    """Calculate median intensity of Red and Green channels."""
    input:
        red=lambda wc: cfg.expand(
            cfg.config.user_files.idat_pattern.red,
            query=f"Sample_ID == '{wc.Sample_ID}'",
        ),
        green=lambda wc: cfg.expand(
            cfg.config.user_files.idat_pattern.green,
            query=f"Sample_ID == '{wc.Sample_ID}'",
        ),
    output:
        temp(
            "sample_level/per_sample_median_idat_intensity/{Sample_ID}.{SentrixBarcode_A}.{SentrixPosition_A}.csv"
        ),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
    group:
        "per_sample_median_idat_intensity"
    conda:
        cfg.conda("illuminaio")
    script:
        "../scripts/median_idat_intensity.R"


rule agg_median_idat_intensity:
    input:
        cfg.expand(rules.per_sample_median_idat_intensity.output[0]),
    output:
        "sample_level/contamination/median_idat_intensity.csv",
    resources:
        mem_gb=lambda wildcards, attempt: attempt * 4,
        tim_hr=lambda wildcards, attempt: attempt ** 2,
    run:
        pd.concat([pd.read_csv(file_name) for file_name in input]).to_csv(output[0], index=False)


# -------------------------------------------------------------------------------
# Convert GTC to ADPC
# -------------------------------------------------------------------------------
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
            cfg.config.user_files.gtc_pattern,
            query=f"Sample_ID == '{wc.Sample_ID}'",
        )[0],
        bpm_file=cfg.config.reference_files.illumina_manifest_file,
    output:
        temp("sample_level/per_sample_adpc/{Sample_ID}.adpc.bin"),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
    group:
        "per_sample_gtc_to_adpc"
    script:
        "../scripts/gtc2adpc.py"


# -------------------------------------------------------------------------------
# Estimate Contamination
# -------------------------------------------------------------------------------
rule per_sample_contamination_verifyIDintensity:
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
        abf=rules.pull_b_allele_freq_from_1kg.output.abf_file,
    params:
        snps=cfg.config.num_snps,
    output:
        temp("sample_level/per_sample_contamination_test/{Sample_ID}.contam.out"),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
    group:
        "per_sample_contamination_verifyIDintensity"
    conda:
        cfg.conda("verifyidintensity")
    shell:
        "verifyIDintensity -m {params.snps} -n 1 -b {input.abf} -v -p -i {input.adpc} > {output}"


rule agg_contamination:
    input:
        cfg.expand(rules.per_sample_contamination_verifyIDintensity.output),
    output:
        "sample_level/contamination/verifyIDintensity.csv",
    resources:
        mem_gb=lambda wildcards, attempt: attempt * 4,
        tim_hr=lambda wildcards, attempt: attempt ** 2,
    run:
        from cgr_gwas_qc.parsers import verifyidintensity

        pd.concat([verifyidintensity.read(sample) for sample in input], ignore_index=True).to_csv(
            output[0], index=False
        )
