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
# PHONY Rules
################################################################################
# NOTE: Because of job grouping it is cleaner to wrap the various CLI utilities in
# their own python script. This makes using conda more complicated. Instead of
# installing all of the python dependencies in each environment, I will just
# pass the conda environment and activate it internal. However, we want to make
# sure that snakemake creates the environments, so these PHONY rules will make
# sure that the conda env exists.
localrules:
    illuminaio_conda,
    verifyidintensity_conda,


rule illuminaio_conda:
    output:
        temp(".illuminaio_env_built"),
    conda:
        cfg.conda("illuminaio")
    shell:
        "touch {output[0]}"


rule verifyidintensity_conda:
    output:
        temp(".verifyidintensity_env_built"),
    conda:
        cfg.conda("verifyidintensity")
    shell:
        "touch {output[0]}"


################################################################################
# Workflow Rules
################################################################################
if config.get("cluster_mode", False):

    rule grouped_median_idat_intensity:
        """Calculate median intensity of Red and Green channels."""
        input:
            sample_sheet_csv="cgr_sample_sheet.csv",
            _=rules.illuminaio_conda.output[0],
        params:
            grp="{grp}",
            idat_red_pattern=lambda _: cfg.config.user_files.idat_pattern.red,
            idat_green_pattern=lambda _: cfg.config.user_files.idat_pattern.green,
            rscript=cfg.scripts("median_idat_intensity.R"),
            conda_env=cfg.conda("illuminaio"),
            notemp=config.get("notemp", False),
        output:
            temp("sample_level/contamination/{grp}/median_idat_intensity.csv"),
        threads: 12
        resources:
            mem_mb=lambda wildcards, attempt: 1024 * 12 * attempt,
            time_hr=lambda wildcards, attempt: 4 * attempt,
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

    rule grouped_contamination:
        """Converts a sample's GTC/BPM to an Illumina ADPC.BIN.

        This is the format required by ``verifyIDintensity``. The script also
        runs some sanity checks (intensities and normalized intensities > 0;
        genotypes are one of {0, 1, 2, 3}) while processing each file.

        .. warning::
            This is a submission hot spot creating 1 job per sample.

        Find contaminated samples using allele intensities.

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
            sample_sheet_csv="cgr_sample_sheet.csv",
            bpm_file=cfg.config.reference_files.illumina_manifest_file,
            abf_file=rules.pull_b_allele_freq_from_1kg.output.abf_file,
            _=rules.verifyidintensity_conda.output[0],
        params:
            grp="{grp}",
            gtc_pattern=lambda _: cfg.config.user_files.gtc_pattern,
            snps=cfg.config.num_snps,
            conda_env=cfg.conda("verifyidintensity"),
            notemp=config.get("notemp", False),
        output:
            temp("sample_level/contamination/{grp}/verifyIDintensity.csv"),
        threads: 12
        resources:
            mem_mb=lambda wildcards, attempt: 1024 * 12 * attempt,
            time_hr=lambda wildcards, attempt: 4 * attempt,
        script:
            "../scripts/grouped_contamination.py"

    rule agg_verifyidintensity:
        input:
            expand(rules.grouped_contamination.output[0], grp=cfg.cluster_groups),
        output:
            "sample_level/contamination/verifyIDintensity.csv",
        resources:
            mem_gb=lambda wildcards, attempt: attempt * 4,
            tim_hr=lambda wildcards, attempt: attempt ** 2,
        script:
            "../scripts/agg_verifyidintensity.py"


else:

    rule per_sample_median_idat_intensity:
        """Calculate median intensity of Red and Green channels."""
        input:
            red=lambda wc: cfg.expand(
                cfg.config.user_files.idat_pattern.red,
                query=f"Sample_ID == '{wc.Sample_ID}'",
            )[0],
            green=lambda wc: cfg.expand(
                cfg.config.user_files.idat_pattern.green,
                query=f"Sample_ID == '{wc.Sample_ID}'",
            )[0],
            _=rules.illuminaio_conda.output[0],
        params:
            sample_id="{Sample_ID}",
            rscript=cfg.scripts("median_idat_intensity.R"),
            conda_env=cfg.conda("illuminaio"),
        output:
            temp("sample_level/contamination/per_sample_median_idat_intensity/{Sample_ID}.csv"),
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
        script:
            "../scripts/median_idat_intensity.py"

    rule agg_median_idat_intensity:
        input:
            cfg.expand(
                rules.per_sample_median_idat_intensity.output[0],
                query="not is_missing_idats",
            ),
        output:
            "sample_level/contamination/median_idat_intensity.csv",
        resources:
            mem_gb=lambda wildcards, attempt: attempt * 4,
            tim_hr=lambda wildcards, attempt: attempt ** 2,
        script:
            "../scripts/agg_median_idat_intensity.py"

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
        script:
            "../scripts/gtc2adpc.py"

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
            _=rules.verifyidintensity_conda.output[0],
        params:
            snps=cfg.config.num_snps,
            conda_env=cfg.conda("verifyidintensity"),
        output:
            temp("sample_level/per_sample_contamination_test/{Sample_ID}.contam.out"),
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
        script:
            "../scripts/verifyidintensity.py"

    rule agg_verifyidintensity:
        input:
            cfg.expand(
                rules.per_sample_contamination_verifyIDintensity.output[0],
                query="not is_missing_gtc",
            ),
        output:
            "sample_level/contamination/verifyIDintensity.csv",
        resources:
            mem_gb=lambda wildcards, attempt: attempt * 4,
            tim_hr=lambda wildcards, attempt: attempt ** 2,
        script:
            "../scripts/agg_verifyidintensity.py"
