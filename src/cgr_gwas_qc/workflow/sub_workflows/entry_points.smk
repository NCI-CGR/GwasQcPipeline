"""Entry points into the QC workflow.

This module contains all the different data entry points into the QC
workflow. The most common case is a set of sample level GTC files provided by
the user.

All entry points should create:

    - sample_level/samples.bed
    - sample_level/samples.bim
    - sample_level/samples.fam
"""

from cgr_gwas_qc import load_config
from cgr_gwas_qc.parsers import sample_sheet
from math import ceil

cfg = load_config()


################################################################################
# Entry Points Targets
################################################################################
targets = [
    "sample_level/samples.bed",
    "sample_level/samples.bim",
    "sample_level/samples.fam",
]


rule all_entry_points:
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
    plink_conda,


rule plink_conda:
    output:
        temp(".plink_env_built"),
    conda:
        cfg.conda("plink2")
    shell:
        "touch {output[0]}"


################################################################################
# Imports
################################################################################
module bcf_module:
    snakefile:
        cfg.modules("bcf")
    config:
        {}


################################################################################
# Workflow Rules
################################################################################
if cfg.config.user_files.gtc_pattern:

    def _get_gtc(wildcards):
        return cfg.expand(
            cfg.config.user_files.gtc_pattern,
            query=f"Sample_ID == '{wildcards.Sample_ID}'",
        )[0]

    ################################################################################
    # GTC To Plink
    ################################################################################

    if not cfg.config.workflow_params.convert_gtc2bcf:
        if config.get("cluster_mode", False):

            rule grouped_gtc_to_bed:
                input:
                    sample_sheet_csv="cgr_sample_sheet.csv",
                    bpm_file=cfg.config.reference_files.illumina_manifest_file,
                    _=rules.plink_conda.output[0],
                params:
                    grp="{grp}",
                    gtc_pattern=lambda wc: cfg.config.user_files.gtc_pattern,
                    strand=cfg.config.software_params.strand,
                    out_prefix="sample_level/{grp}/samples",
                    conda_env=cfg.conda("plink2"),
                    notemp=config.get("notemp", False),
                output:
                    bed=temp("sample_level/{grp}/samples.bed"),
                    bim=temp("sample_level/{grp}/samples.bim"),
                    fam=temp("sample_level/{grp}/samples.fam"),
                    nosex=temp("sample_level/{grp}/samples.nosex"),
                threads: 12
                resources:
                    mem_mb=lambda wildcards, attempt: 1024 * 12 * attempt,
                    time_hr=lambda wildcards, attempt: 4 * attempt,
                benchmark:
                    "benchmarks/gtc_to_bed.{grp}" + ".tsv"
                script:
                    "../scripts/grouped_gtc_to_bed.py"

            rule merge_grouped_beds:
                input:
                    bed=expand("sample_level/{grp}/samples.bed", grp=cfg.cluster_groups),
                    bim=expand("sample_level/{grp}/samples.bim", grp=cfg.cluster_groups),
                    fam=expand("sample_level/{grp}/samples.fam", grp=cfg.cluster_groups),
                    _=rules.plink_conda.output[0],
                params:
                    out_prefix="sample_level/samples",
                    conda_env=cfg.conda("plink2"),
                    notemp=config.get("notemp", False),
                output:
                    bed="sample_level/samples.bed",
                    bim="sample_level/samples.bim",
                    fam="sample_level/samples.fam",
                    nosex="sample_level/samples.nosex",
                log:
                    "sample_level/samples.log",
                threads: 8
                resources:
                    mem_mb=lambda wildcards, attempt: 1024 * 8 * attempt,
                    time_hr=lambda wildcards, attempt: 4 * attempt,
                benchmark:
                    "benchmarks/plink_merge" + ".tsv"
                script:
                    "../scripts/plink_merge.py"

        else:

            rule per_sample_gtc_to_ped:
                """Converts a sample's GTC/BPM to PED/MAP.

                .. warning::
                    This is a submission hot spot creating 1 job per sample.
                """
                input:
                    gtc=_get_gtc,
                    bpm=cfg.config.reference_files.illumina_manifest_file,
                params:
                    strand=cfg.config.software_params.strand,
                output:
                    ped=temp("sample_level/per_sample_plink_files/{Sample_ID}.ped"),
                    map_=temp("sample_level/per_sample_plink_files/{Sample_ID}.map"),
                resources:
                    mem_mb=lambda wildcards, attempt: attempt * 1024,
                script:
                    "../scripts/gtc2plink.py"

            rule merge_gtc_sample_peds:
                """Merges multiple samples using ``plink --merge-list``.

                Merge and convert sample(s) into an aggregated binary ``plink``
                format.
                """
                input:
                    ped=cfg.expand(rules.per_sample_gtc_to_ped.output.ped),
                    map_=cfg.expand(rules.per_sample_gtc_to_ped.output.map_),
                    _=rules.plink_conda.output[0],
                params:
                    out_prefix="sample_level/samples",
                    conda_env=cfg.conda("plink2"),
                    notemp=config.get("notemp", False),
                output:
                    bed="sample_level/samples.bed",
                    bim="sample_level/samples.bim",
                    fam="sample_level/samples.fam",
                    nosex="sample_level/samples.nosex",
                log:
                    "sample_level/samples.log",
                threads: 8
                resources:
                    mem_mb=lambda wildcards, attempt: 1024 * 8 * attempt,
                    time_hr=lambda wildcards, attempt: 4 * attempt,
                script:
                    "../scripts/plink_merge.py"

    ################################################################################
    # GTC To BCF
    ################################################################################

    if cfg.config.workflow_params.convert_gtc2bcf:

        localrules:
            gtc2bcf_conda,
            write_gtc_pathlist,

        rule gtc2bcf_conda:
            output:
                temp(".bcftools-gtc2vcf-plugin_env_built"),
            conda:
                cfg.conda("bcftools-gtc2vcf-plugin")
            shell:
                "touch {output[0]}"

        if config.get("cluster_mode", False) and len(cfg.cluster_groups) > 2:

            def _get_n_samples(wildcards):
                return len(cfg.ss.query(f'cluster_group=="{wildcards.grp}"'))

            use rule write_gtc_pathlist from bcf_module with:
                params:
                    grp=cfg.cluster_groups,
                output:
                    temp("sample_level/{grp}/gtc.tsv"),

            use rule gtc_to_bcf from bcf_module with:
                input:
                    gtcs=rules.write_gtc_pathlist.output[0],
                benchmark:
                    "benchmarks/gtc_to_bcf.{grp}" + ".tsv"
                resources:
                    time_hr=lambda wc: ceil(
                        ((_get_n_samples(wc) + 1) * (cfg.config.num_snps * 3e-6)) / 3600
                    )
                    + 1,
                    mem_mb=lambda wc: ceil(
                        (_get_n_samples(wc) * (cfg.config.num_snps * 1.06e-6))
                        + (cfg.config.num_snps * 2e-3)
                    )
                    + 200,
                output:
                    bcf="sample_level/{grp}/samples.bcf",

            rule merge_gtc_to_bcf_batches:
                """Merges grouped bcf files to single merged bcf"""
                input:
                    bcf_batches=expand("sample_level/{grp}/samples.bcf", grp=cfg.cluster_groups),
                output:
                    bcf="sample_level/samples.bcf",
                conda:
                    cfg.conda("bcftools-gtc2vcf-plugin")
                benchmark:
                    (
                        "benchmarks/merge_gtc_to_bcf_batches"
                        + "."
                        + str(len(cfg.cluster_groups))
                        + "."
                        + str(len(cfg.ss))
                        + ".tsv"
                    )
                threads: 44
                resources:
                    time_hr=ceil((len(cfg.ss) * 0.2) / 3600 + 1),
                    mem_mb=len(cfg.cluster_groups) * 70,
                shell:
                    "bcftools merge --threads {threads} --merge none {input.bcf_batches} -Ob -o {output}"

            use rule convert_bcf_to_plink_bed from bcf_module with:
                input:
                    bcf=rules.merge_gtc_to_bcf_batches.output[0],

        else:

            use rule write_gtc_pathlist from bcf_module with:
                params:
                    grp="",
                output:
                    "sample_level/gtc.tsv",

            use rule gtc_to_bcf from bcf_module with:
                output:
                    bcf="sample_level/samples.bcf",

            use rule convert_bcf_to_plink_bed from bcf_module with:
                input:
                    bcf=rules.gtc_to_bcf.output[0],

elif cfg.config.user_files.ped and cfg.config.user_files.map:

    ################################################################################
    # Aggregated PED/MAPs
    ################################################################################
    rule plink_make_bed:
        """Converts aggregated PED/MAP into binary ``plink`` format."""
        input:
            ped=cfg.config.user_files.ped,
            map_=cfg.config.user_files.map,
        params:
            prefix="sample_level/samples",
        output:
            bed="sample_level/samples.bed",
            bim="sample_level/samples.bim",
            fam="sample_level/samples.fam",
            nosex="sample_level/samples.nosex",
        log:
            "sample_level/samples.log",
        conda:
            cfg.conda("plink2")
        resources:
            mem_mb=lambda wildcards, attempt: 1024 * 8 * attempt,
        shell:
            "plink "
            "--ped {input.ped} "
            "--map {input.map_} "
            "--make-bed "
            "--out {params.prefix} "
            "--memory {resources.mem_mb}"

elif cfg.config.user_files.bed and cfg.config.user_files.bim and cfg.config.user_files.fam:

    ################################################################################
    # Aggregated BED/BIM/FAM
    ################################################################################
    localrules:
        link_binary_plink_to_expected_location,

    rule link_binary_plink_to_expected_location:
        """Creates symbolic links BED/BIM/FAM.

        BED/BIM/FAM files are expected to be at
        ``plink_start/samples.{bed,bim,fam}``. Since the user is providing
        these files, we are using symbolic links to place them in the correct
        location.
        """
        input:
            bed=cfg.config.user_files.bed,
            bim=cfg.config.user_files.bim,
            fam=cfg.config.user_files.fam,
        output:
            bed="sample_level/samples.bed",
            bim="sample_level/samples.bim",
            fam="sample_level/samples.fam",
        shell:
            "ln -s $(readlink -f {input.bed}) {output.bed} && "
            "ln -s $(readlink -f {input.bim}) {output.bim} && "
            "ln -s $(readlink -f {input.fam}) {output.fam}"

elif cfg.config.user_files.bcf:

    ################################################################################
    # Aggregated BCF
    ################################################################################
    localrules:
        convert_bcf_to_plink_bed,

    use rule convert_bcf_to_plink_bed from bcf_module with:
        input:
            bcf=cfg.config.user_files.bcf,
