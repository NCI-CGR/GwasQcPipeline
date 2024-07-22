import pandas as pd

from cgr_gwas_qc import load_config

cfg = load_config()


################################################################################
# Intensity Targets
################################################################################
targets = [
    "sample_level/contamination/median_idat_intensity.csv",
]


rule all_contamination:
    input:
        targets,


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


rule illuminaio_conda:
    output:
        temp(".illuminaio_env_built"),
    conda:
        cfg.conda("illuminaio")
    shell:
        "touch {output[0]}"


################################################################################
# Workflow Rules
################################################################################
if cfg.config.user_files.bcf:
    if config.get("cluster_mode", False):

        localrules:
            agg_median_idat_intensity,

        rule grouped_median_intensity_from_vcf:
            """Calculate median intensity from raw intensities using VCF/BCF input."""
            input:
                sample_sheet_csv="cgr_sample_sheet.csv",
                vcf_file=cfg.config.user_files.bcf,
            params:
                grp="{grp}",
                notemp=config.get("notemp", False),
            output:
                temp("sample_level/contamination/{grp}/median_idat_intensity.csv"),
            threads: 12
            resources:
                mem_mb=lambda wildcards, attempt: 1024 * 12 * attempt,
                time_hr=lambda wildcards, attempt: 4 * attempt,
            script:
                "../scripts/grouped_intensity_from_vcf.py"

        rule agg_median_idat_intensity:
            input:
                expand(rules.grouped_median_intensity_from_vcf.output[0], grp=cfg.cluster_groups),
            params:
                notemp=config.get("notemp", False),
            output:
                "sample_level/contamination/median_idat_intensity.csv",
            resources:
                mem_gb=lambda wildcards, attempt: attempt * 4,
                time_hr=lambda wildcards, attempt: attempt**2,
            script:
                "../scripts/agg_median_idat_intensity.py"

    else:

        rule per_sample_median_intensity_from_vcf:
            """Calculate median intensity from raw intensities using VCF/BCF input."""
            input:
                vcf_file=cfg.config.user_files.bcf,
            params:
                sample_id="{Sample_ID}",
            output:
                temp("sample_level/contamination/per_sample_median_idat_intensity/{Sample_ID}.csv"),
            resources:
                mem_mb=lambda wildcards, attempt: attempt * 1024,
            script:
                "../scripts/median_intensity_from_vcf.py"

        rule agg_median_idat_intensity:
            input:
                cfg.expand(rules.per_sample_median_intensity_from_vcf.output[0]),
            output:
                "sample_level/contamination/median_idat_intensity.csv",
            resources:
                mem_gb=lambda wildcards, attempt: attempt * 4,
                time_hr=lambda wildcards, attempt: attempt**2,
            script:
                "../scripts/agg_median_idat_intensity.py"

else:
    if config.get("cluster_mode", False):

        localrules:
            agg_median_idat_intensity,

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
                time_hr=lambda wildcards, attempt: attempt**2,
            script:
                "../scripts/agg_median_idat_intensity.py"

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
                cfg.expand(rules.per_sample_median_idat_intensity.output[0]),
            output:
                "sample_level/contamination/median_idat_intensity.csv",
            resources:
                mem_gb=lambda wildcards, attempt: attempt * 4,
                time_hr=lambda wildcards, attempt: attempt**2,
            script:
                "../scripts/agg_median_idat_intensity.py"
