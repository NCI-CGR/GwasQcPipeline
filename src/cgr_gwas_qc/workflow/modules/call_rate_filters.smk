
rule call_rate_filter_1:
    """Outputs ``plink`` dataset removing samples and variants below call rate filter 1.

    Call rate filter 1 is usually (1 - 0.8) missing for samples and variants.

    .. resources::
        :memory: 10GB

    """
    input:
        bed="plink_start/samples.bed",
        bim="plink_start/samples.bim",
        fam="plink_start/samples.fam",
    params:
        mind=1 - cfg.software_params.samp_cr_1,
        geno=1 - cfg.software_params.snp_cr_1,
    output:
        bed="plink_filter_call_rate_1/samples.bed",
        bim="plink_filter_call_rate_1/samples.bim",
        fam="plink_filter_call_rate_1/samples.fam",
    log: "plink_filter_call_rate_1/samples.log",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    threads: 20
    resources:
        mem=10000,
    shell:
        "plink "
        "--bfile $(dirname {input[0]}) "
        "--mind {params.mind} "
        "--geno {params.geno} "
        "--make-bed "
        "--memory {resources.mem} "
        "--out $(dirname {output[0]})"


rule call_rate_filter_2:
    """Outputs ``plink`` dataset removing samples and variants below call rate filter 2.

    Call rate filter 2 is usually (1 - 0.95) missing for samples and variants.

    .. resources::
        :memory: 10GB

    """
    input:
        bed="plink_filter_call_rate_1/samples.bed",
        bim="plink_filter_call_rate_1/samples.bim",
        fam="plink_filter_call_rate_1/samples.fam",
    params:
        mind=1 - cfg.software_params.samp_cr_2,
        geno=1 - cfg.software_params.snp_cr_2,
    output:
        bed="plink_filter_call_rate_2/samples.bed",
        bim="plink_filter_call_rate_2/samples.bim",
        fam="plink_filter_call_rate_2/samples.fam",
    log: "plink_filter_call_rate_2/samples.log",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    threads: 20
    resources:
        mem=10000,
    shell:
        "plink "
        "--bfile $(dirname {input[0]}) "
        "--mind {params.mind} "
        "--geno {params.geno} "
        "--make-bed "
        "--memory {resources.mem} "
        "--out $(dirname {output[0]})"
