import pandas as pd

from cgr_gwas_qc import load_config

cfg = load_config()


################################################################################
# Contamination Targets
################################################################################
targets = [
    "sample_level/contamination/median_idat_intensity.csv",
    "sample_level/contamination/verifyIDintensity.csv",
]


rule all_contamination:
    input:
        targets,


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
if config.get("cluster_mode", False):

    rule illuminaio_conda:
        """Fake rule to ensure the illuminaio conda env is built.

        The grouped median idat intensity requires this conda env, but uses it
        internally. I want to make sure that snakemake builds the conda env.
        """
        output:
            temp(".illuminaio_env_built"),
        conda:
            cfg.conda("illuminaio")
        shell:
            "touch {output[0]}"

    rule grouped_median_idat_intensity:
        """Calculate median intensity of Red and Green channels."""
        input:
            sample_sheet_csv="cgr_sample_sheet.csv",
            _=rules.illuminaio_conda.output[0],
        params:
            grp="{grp}",
            idat_red_pattern=lambda x: cfg.config.user_files.idat_pattern.red,
            idat_green_pattern=lambda x: cfg.config.user_files.idat_pattern.green,
            conda_env=cfg.conda("illuminaio"),
            notemp=config.get("notemp", False),
        output:
            temp("sample_level/contamination/{grp}/median_idat_intensity.csv"),
        threads: 12
        resources:
            mem_mb=lambda wildcards, attempt: 1024 * 12 * attempt,
            time_hr=lambda wildcards, attempt: 4 * attempt,
        group:
            "per_sample_median_idat_intensity"
        script:
            "../scripts/grouped_median_idat_intensity.py"

    rule agg_median_idat_intensity:
        input:
            expand(rules.grouped_median_idat_intensity.output[0], grp=cfg.cluster_groups),
        params:
            notemp=config.get("notemp", False),
        output:
            "sample_level/contamination/median_idat_intensity.csv",
        resources:
            mem_gb=lambda wildcards, attempt: attempt * 4,
            tim_hr=lambda wildcards, attempt: attempt ** 2,
        script:
            "../scripts/agg_median_idat_intensity.py"


else:

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
        script:
            "../scripts/agg_median_idat_intensity.py"


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
