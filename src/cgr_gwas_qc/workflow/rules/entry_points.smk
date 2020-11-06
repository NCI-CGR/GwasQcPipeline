"""Entry points into the QC workflow.

This module contains all the different data entry points into the QC
workflow. The most common case is a set of sample level GTC files provided by
the user.

All entry points should create:

    - plink_start/samples.bed
    - plink_start/samples.bim
    - plink_start/samples.fam
"""

if cfg.config.user_files.gtc_pattern:
    ################################################################################
    # GTC To Plink
    ################################################################################
    rule gtc_to_ped:
        """Converts a sample's GTC/BPM to PED/MAP.

        .. warning::
            This is a submission hot spot creating 1 job per sample.
        """
        input:
            gtc=lambda wc: cfg.expand(
                cfg.config.user_files.gtc_pattern, query=f"Sample_ID == '{wc.Sample_ID}'"
            )[0],
            bpm=cfg.config.reference_files.illumina_manifest_file,
        params:
            strand=cfg.config.software_params.strand,
        output:
            ped=temp("per_sample_plink_files/{Sample_ID}.ped"),
            map_=temp("per_sample_plink_files/{Sample_ID}.map"),
        script:
            "../scripts/gtc2plink.py"

    rule create_gtc_sample_merge_list:
        """Makes a list of samples to merge together.

        Plink needs a list of samples (PEDs, MAPs) to merge together. This is
        a simple space delimited file where each sample is on it's own row.
        """
        input:
            ped=cfg.expand("per_sample_plink_files/{Sample_ID}.ped"),
            map_=cfg.expand("per_sample_plink_files/{Sample_ID}.map"),
        output:
            "plink_start/mergeList.txt",
        group:
            "merge_gtc_sample_peds"
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
            ped=cfg.expand("per_sample_plink_files/{Sample_ID}.ped"),
            map_=cfg.expand("per_sample_plink_files/{Sample_ID}.map"),
            merge_list=rules.create_gtc_sample_merge_list.output[0],
        params:
            prefix="plink_start/samples",
        output:
            bed="plink_start/samples.bed",
            bim="plink_start/samples.bim",
            fam="plink_start/samples.fam",
            nosex="plink_start/samples.nosex",
        log: "plink_start/samples.log",
        group:
            "merge_gtc_sample_peds"
        envmodules:
            cfg.envmodules("plink2"),
        conda:
            "../conda/plink2.yml"
        threads: 20
        resources:
            mem=24000,
        shell:
            "plink "
            "--merge-list {input.merge_list} "
            "--make-bed "
            "--out {params.prefix} "
            "--threads {threads} "
            "--memory {resources.mem}"

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
            prefix="plink_start/samples",
        output:
            bed="plink_start/samples.bed",
            bim="plink_start/samples.bim",
            fam="plink_start/samples.fam",
            nosex="plink_start/samples.nosex",
        envmodules:
            cfg.envmodules("plink2"),
        conda:
            "../conda/plink2.yml"
        resources:
            mem=10000,
        shell:
            "plink "
            "--ped {input.ped} "
            "--map {input.map_} "
            "--make-bed "
            "--out {params.prefix} "
            "--memory {resources.mem}"

elif cfg.config.user_files.bed and cfg.config.user_files.bim and cfg.config.user_files.fam:
    ################################################################################
    # Aggregated BED/BIM/FAM
    ################################################################################
    rule link_binary_plink_to_expected_loction:
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
            bed="plink_start/samples.bed",
            bim="plink_start/samples.bim",
            fam="plink_start/samples.fam",
        shell:
            "ln -s $(readlink -f {input.bed}) {output.bed} && "
            "ln -s $(readlink -f {input.bim}) {output.bim} && "
            "ln -s $(readlink -f {input.fam}) {output.fam}"
