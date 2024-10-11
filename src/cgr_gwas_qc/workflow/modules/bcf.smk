from cgr_gwas_qc import load_config
from math import ceil

cfg = load_config()

############################
# General BCF import rules
###############################


def _create_unknown_sex(wildcards):
    path_to_unknown_sex_lst = temp("unknown_sex.lst")
    cfg.ss.rename(columns={"Sample_ID": "#IID"}).assign(Sex="U").to_csv(
        path_to_unknown_sex_lst, sep="\t", index=False
    )
    return path_to_unknown_sex_lst


rule convert_bcf_to_plink_bed:
    """Converts BCF to plink BED file

    Path to aggregated BCF file is expected in user_files in config. The expected
    BCF should be created using BCFtools/gtc2vcf. The BCF will be converted to BED
    file for processing.
    """
    input:
        bcf=cfg.config.user_files.bcf,
    params:
        prefix="sample_level/samples",
        unknown_sex=_create_unknown_sex,
    output:
        bed="sample_level/samples.bed",
        bim="sample_level/samples.bim",
        fam="sample_level/samples.fam",
    log:
        "sample_level/samples.log",
    conda:
        cfg.conda("plink2-0")
    benchmark:
        "benchmarks/convert_bcf_to_plink_bed." + str(len(cfg.ss)) + ".tsv"
    threads: 44
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * 8 * attempt,
        time_hr=ceil((0.11 * len(cfg.ss)) / 3600),
    shell:
        "plink2 --allow-extra-chr 0 --keep-allele-order --double-id --bcf {input.bcf} --update-sex {params.unknown_sex} --output-chr 26 --split-par hg38 --make-pgen --out sample_level/bcf2plink  --memory {resources.mem_mb} --threads {threads} ;"
        "plink2 --pfile sample_level/bcf2plink --make-pgen --sort-vars --out sample_level/bcf2plink-sorted --threads {threads} --memory {resources.mem_mb}  ;"
        "plink2 --pfile sample_level/bcf2plink-sorted --make-bed --out sample_level/samples --threads {threads} --memory {resources.mem_mb} ;"
        "rm sample_level/bcf2plink.{{pgen,psam,pvar,log}} sample_level/bcf2plink-sorted.{{pgen,psam,pvar,log}}"


rule write_gtc_pathlist:
    """Writes a TSV for GTC path list as input to be used in gtc_to_bcf rule
    """
    params:
        grp="",
    output:
        "sample_level/gtc.tsv",
    run:
        if params.grp == "":
            gtcList = cfg.expand(cfg.config.user_files.gtc_pattern)
        else:
            gtcList = cfg.expand(
                cfg.config.user_files.gtc_pattern, query=f'cluster_group=="{wildcards.grp}"'
            )
        with open(output[0], "w") as f:
            for line in gtcList:
                f.write(f"{line}\n")


def _get_add_params_for_gtc2bcf(wildcards):
    if cfg.config.reference_files.illumina_csv_bpm:
        return (
            "--csv "
            + str(cfg.config.reference_files.illumina_csv_bpm)
            + " "
            + str(cfg.config.workflow_params.additional_params_for_gtc2bcf)
        )
    else:
        return str(cfg.config.workflow_params.additional_params_for_gtc2bcf)


rule gtc_to_bcf:
    """Converts a list of GTC files to BCF.
    """
    input:
        gtcs=rules.write_gtc_pathlist.output[0],
    params:
        additional_params=_get_add_params_for_gtc2bcf,
        bpm=cfg.config.reference_files.illumina_manifest_file,
        reference_fasta=cfg.config.reference_files.reference_fasta,
    output:
        bcf="sample_level/samples.bcf",
    threads: 4
    benchmark:
        "benchmarks/gtc_to_bcf." + str(len(cfg.ss)) + ".tsv"
    conda:
        cfg.conda("bcftools-gtc2vcf-plugin")
    resources:
        time_hr=ceil(((len(cfg.ss) + 1) * (cfg.config.num_snps * 3e-6)) / 3600) + 1,
        mem_mb=ceil((len(cfg.ss) * (cfg.config.num_snps * 1.06e-6)) + (cfg.config.num_snps * 2e-3))
        + 200,
    shell:
        """
        bcftools +gtc2vcf --threads {threads} --gtcs {input.gtcs} --bpm {params.bpm} --fasta-ref {params.reference_fasta} {params.additional_params} -Ou | bcftools sort -Ou -T ./bcftools. | bcftools norm --no-version -Ob --check-ref x -f {params.reference_fasta} --multiallelics -any --write-index --output {output.bcf}
        """
