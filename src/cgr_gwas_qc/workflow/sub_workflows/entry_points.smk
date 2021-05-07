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

cfg = load_config()


################################################################################
# All Targets
################################################################################
rule all_entry_points:
    input:
        "sample_level/samples.bed",
        "sample_level/samples.bim",
        "sample_level/samples.fam",


################################################################################
# Workflow Rules
################################################################################
if cfg.config.user_files.gtc_pattern:

    ################################################################################
    # GTC To Plink
    ################################################################################
    rule per_sample_gtc_to_ped:
        """Converts a sample's GTC/BPM to PED/MAP.

        .. warning::
            This is a submission hot spot creating 1 job per sample.
        """
        input:
            gtc=lambda wc: cfg.expand(
                cfg.config.user_files.gtc_pattern, query=f"Sample_ID == '{wc.Sample_ID}'",
            )[0],
            bpm=cfg.config.reference_files.illumina_manifest_file,
        params:
            strand=cfg.config.software_params.strand,
        output:
            ped=temp("sample_level/per_sample_plink_files/{Sample_ID}.ped"),
            map_=temp("sample_level/per_sample_plink_files/{Sample_ID}.map"),
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
        group:
            "per_sample_gtc_to_ped"
        script:
            "../scripts/gtc2plink.py"

    rule create_gtc_sample_merge_list:
        """Makes a list of samples to merge together.

        Plink needs a list of samples (PEDs, MAPs) to merge together. This is
        a simple space delimited file where each sample is on it's own row.
        """
        input:
            ped=cfg.expand(rules.per_sample_gtc_to_ped.output.ped),
            map_=cfg.expand(rules.per_sample_gtc_to_ped.output.map_),
        output:
            temp("sample_level/initial_mergeList.txt"),
        run:
            with open(output[0], "w") as fh:
                for ped, map_ in zip(input.ped, input.map_):
                    fh.write(f"{ped} {map_}\n")

    rule merge_gtc_sample_peds:
        """Merges multiple samples using ``plink --merge-list``.

        Merge and convert sample(s) into an aggregated binary ``plink``
        format.
        """
        input:
            ped=cfg.expand(rules.per_sample_gtc_to_ped.output.ped),
            map_=cfg.expand(rules.per_sample_gtc_to_ped.output.map_),
            merge_list=rules.create_gtc_sample_merge_list.output[0],
        params:
            prefix="sample_level/samples",
        output:
            bed="sample_level/samples.bed",
            bim="sample_level/samples.bim",
            fam="sample_level/samples.fam",
            nosex="sample_level/samples.nosex",
        log:
            "sample_level/samples.log",
        envmodules:
            cfg.envmodules("plink2"),
        conda:
            cfg.conda("plink2")
        threads: 8
        resources:
            mem_mb=lambda wildcards, attempt: 1024 * 8 * attempt,
            time_hr=lambda wildcards, attempt: 4 * attempt,
        shell:
            "plink "
            "--merge-list {input.merge_list} "
            "--make-bed "
            "--out {params.prefix} "
            "--threads {threads} "
            "--memory {resources.mem_mb}"


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
        envmodules:
            cfg.envmodules("plink2"),
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
