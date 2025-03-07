from cgr_gwas_qc import load_config
from math import ceil

cfg = load_config()

############################
# General BCF import rules
###############################


def _create_unknown_sex(wildcards):
    import tempfile

    path_to_unknown_sex_lst = tempfile.NamedTemporaryFile(delete=False)
    cfg.ss.rename(columns={"Sample_ID": "#IID"}).assign(Sex="U").to_csv(
        path_to_unknown_sex_lst, sep="\t", index=False
    )
    return path_to_unknown_sex_lst.name


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
    threads: workflow.cores
    resources:
        mem_mb=ceil((0.07 * len(cfg.ss))) + 1024,
        time_hr=ceil((0.22 * len(cfg.ss)) + (9e-4 * cfg.config.num_snps) / 3600),
    shell:
        "plink2 --allow-extra-chr 0 --keep-allele-order --double-id --bcf {input.bcf} --vcf-filter --update-sex {params.unknown_sex} --output-chr 26 --split-par hg38 --make-pgen --out sample_level/bcf2plink  --memory {resources.mem_mb} --threads {threads} ;"
        "plink2 --pfile sample_level/bcf2plink --make-pgen --sort-vars --out sample_level/bcf2plink-sorted --threads {threads} --memory {resources.mem_mb}  ;"
        "plink2 --pfile sample_level/bcf2plink-sorted --make-bed --out sample_level/samples --threads {threads} --memory {resources.mem_mb} ;"
        "rm sample_level/bcf2plink.{{pgen,psam,pvar,log}} sample_level/bcf2plink-sorted.{{pgen,psam,pvar,log}} {params.unknown_sex}"


rule write_gtc_pathlist:
    """Writes a TSV for GTC path list as input to be used in gtc_to_bcf rule
    """
    params:
        grp="",
        pattern=lambda wc: cfg.config.user_files.gtc_pattern,
    output:
        "sample_level/gtc.tsv",
    run:
        if params.grp == "":
            gtcList = cfg.expand(params.pattern, query="is_missing_gtc==False")
        else:
            params.pattern = expand(params.pattern, grp=wildcards.grp, allow_missing=True)
            gtcList = cfg.expand(
                params.pattern,
                query='cluster_group=="{grp}"&is_missing_gtc==False'.format(grp=wildcards.grp),
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
        gtcs="sample_level/gtc.tsv",
    params:
        additional_params=_get_add_params_for_gtc2bcf,
        bpm=cfg.config.reference_files.illumina_manifest_file,
        reference_fasta=cfg.config.reference_files.reference_fasta,
        gtc2vcf_location=cfg.SRC_DIR.as_posix() + "/parsers/bcftools-plugins/gtc2vcf.so",
    output:
        bcf="sample_level/samples.bcf",
    threads: 4
    benchmark:
        "benchmarks/gtc_to_bcf." + str(len(cfg.ss)) + ".tsv"
    conda:
        cfg.conda("bcftools")
    resources:
        time_hr=ceil(((len(cfg.ss) + 1) * (cfg.config.num_snps * 3e-6)) / 3600) + 1,
        mem_mb=ceil((len(cfg.ss) * (cfg.config.num_snps * 1.06e-6)) + (cfg.config.num_snps * 2e-3))
        + 200,
    shell:
        """
        bcftools +{params.gtc2vcf_location} --threads {threads} --gtcs {input.gtcs} --bpm {params.bpm} --fasta-ref {params.reference_fasta} {params.additional_params} -Ou | bcftools sort -Ou -T ./bcftools. | bcftools norm --no-version -Ou --check-ref x -f {params.reference_fasta} --multiallelics -any |
        bcftools filter --exclude 'REF==ALT|INFO/INTENSITY_ONLY=1' --soft-filter 'int_only' -Ob --write-index --output {output.bcf}
        """


rule symlink_bcf:
    """Symlinks BCF file to user specified location
    """
    input:
        bcf=cfg.config.user_files.bcf,
    output:
        "sample_level/samples.bcf",
    shell:
        "ln -s {input.bcf} {output}"
