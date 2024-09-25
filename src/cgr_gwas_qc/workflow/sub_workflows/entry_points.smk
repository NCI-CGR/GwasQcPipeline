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
# Workflow Rules
################################################################################
if cfg.config.user_files.gtc_pattern:
    ################################################################################
    # GTC To Plink
    ################################################################################
    def _get_gtc():
        sample_sheet_csv = "cgr_sample_sheet.csv"
        ss = sample_sheet.read(sample_sheet_csv)
        gtcs = []
        for _, record in ss.iterrows():
            gtcs.append(cfg.config.user_files.gtc_pattern.format(**record.to_dict()))
        with open("gtcs.tsv", "w") as file:
            for gtc in gtcs:
                file.write("%s\n" % gtc)

    _get_gtc()

    rule gtc_to_vcf:
        input:
            sample_sheet_csv="cgr_sample_sheet.csv",
            gtcs="gtcs.tsv",
            bpm=cfg.config.reference_files.illumina_manifest_file,
            reference_fasta=cfg.config.reference_files.reference_fasta,
        output:
            vcf=temp("sample_level/samples.vcf"),
        threads: 12
        resources:
            mem_mb=lambda wildcards, attempt: 1024 * 12 * attempt,
            time_hr=lambda wildcards, attempt: 120,
        conda:
            cfg.conda("bcftools-gtc2vcf-plugin")
        shell:
            "bcftools +gtc2vcf --gtcs {input.gtcs} --bpm {input.bpm} --fasta-ref {input.reference_fasta} --output {output.vcf} --use-gtc-sample-names"

    rule filter_missing_allele_snps:
        input:
            vcf=rules.gtc_to_vcf.output.vcf,
        output:
            vcf=temp("sample_level/samples_filtered.vcf"),
        shell:
            "grep -vP '\t\.\t\.\t\.' {input.vcf} > {output.vcf}"

    rule vcf_to_bed:
        input:
            vcf=rules.filter_missing_allele_snps.output.vcf,
        params:
            out_prefix="sample_level/samples",
        output:
            bed="sample_level/samples.bed",
            bim="sample_level/samples.bim",
            fam="sample_level/samples.fam",
            nosex="sample_level/samples.nosex",
        log:
            "sample_level/samples.log",
        conda:
            cfg.conda("plink2")
        threads: 8
        resources:
            mem_mb=lambda wildcards, attempt: 1024 * 8 * attempt,
            time_hr=lambda wildcards, attempt: 4 * attempt,
        shell:
            "plink --vcf {input.vcf} --make-bed --out {params.out_prefix} --double-id"

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

    rule convert_bcf_to_plink_bed:
        """Converts BCF to plink BED file

        Path to aggregated BCF file is expected in user_files in config. The expected
        VCF should be created using BCFtools/gtc2vcf. The BCF will be converted to BED
        file for processing.
        """
        input:
            bcf=cfg.config.user_files.bcf,
        params:
            prefix="sample_level/samples",
        output:
            bed="sample_level/samples.bed",
            bim="sample_level/samples.bim",
            fam="sample_level/samples.fam",
        log:
            "sample_level/samples.log",
        conda:
            cfg.conda("plink2")
        resources:
            mem_mb=lambda wildcards, attempt: 1024 * 8 * attempt,
        shell:
            "plink "
            "--allow-extra-chr 0 --keep-allele-order --double-id "
            "--bcf {input.bcf} "
            "--make-bed "
            "--out {params.prefix} "
            "--memory {resources.mem_mb}"
