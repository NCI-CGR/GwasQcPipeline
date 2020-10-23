
rule plink_sample_qc_stats:
    """Runs ``plink`` statistics.

    ``plink`` will calculate a number of statistics by passing them as a flag.

    .. rubric:: Statistics

        :imiss`: sample-based missing data report
        :lmiss`: variant-based missing data report
        :sexcheck`: Runs X chromosome based sanity checks to verify sex
        :frq`: The minor allele frequency (MAF)
        :hwe`: Checks variants are in HWE

    .. resources::
        :memory: 10GB

    """
    input:
        bed="{prefix}/samples.bed",
        bim="{prefix}/samples.bim",
        fam="{prefix}/samples.fam",
    output:
        imiss="{prefix}/samples{filter}.imiss",
        lmiss="{prefix}/samples{filter}.lmiss",
        sexcheck="{prefix}/samples{filter}.sexcheck",
        frq="{prefix}/samples{filter}.frq",
        hwe="{prefix}/samples{filter}.hwe",
    resources:
        mem=10000,
    shell:
        "plink "
        "--bfile $(dirname {input[0]}) "
        "--freq "
        "--missing "
        "--hardy "
        "--het "
        "--check-sex "
        "--memory {resources.mem} "
        "--out $(dirname {output[0]}"
