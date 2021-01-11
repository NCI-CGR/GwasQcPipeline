""" Sample Concordance: Validate known replicates and identify unknown highly similar samples.

Starts with LD filtering to remove markers that are highly correlated for
better IBS/IBD estimates. Then uses sequence similarity between all pairs of
samples to get a measure of concordance. Specifically, uses the proportion of
shared homozygous sites.
"""


rule plink_ibd:
    """Calculates Identity-by-descent.

    This step is trying to find relationships among the samples by
    calculating IBS/IBD. This calculation is no LD-aware so we are doing the
    LD pruning before. This calculation excludes non-autosomes.
    """
    input:
        bed="snp_filters/ld_prune/samples.bed",
        bim="snp_filters/ld_prune/samples.bim",
        fam="snp_filters/ld_prune/samples.fam",
    params:
        in_prefix="snp_filters/ld_prune/samples",
        out_prefix="sample_concordance/ibd/samples",
        ibd_min=cfg.config.software_params.ibd_pi_hat_min,
        ibd_max=cfg.config.software_params.ibd_pi_hat_max,
    output:
        genome="sample_concordance/ibd/samples.genome",
        nosex="sample_concordance/ibd/samples.nosex",
    log:
        "sample_concordance/ibd/samples.log",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    threads: 20
    resources:
        mem=24000,
    shell:
        "plink "
        "--bfile {params.in_prefix} "
        "--genome full "
        "--min {params.ibd_min} "
        "--max {params.ibd_max} "
        "--threads {threads} "
        "--memory {resources.mem} "
        "--out {params.out_prefix}"


rule sample_concordance:
    """Summarize sample concordance using IBS/IBD.

    Calculates the proportion of shared homozygous markers (IBS2 / (IBS0 + IBS1 + IBS2)) as a
    measure of sample concordance. Then outputs concordance measures for samples that are known
    to be replicates and samples that are thought to be unrelated/independent with a concordance
    > dup_concordance_cutoff (currently 0.95).
    """
    input:
        sample_sheet=cfg.config.sample_sheet,
        imiss="plink_filter_call_rate_2/samples.imiss",
        ibd="sample_concordance/ibd/samples.genome",
    params:
        subject_id_override=cfg.config.workflow_params.subject_id_to_use,
        concordance_threshold=cfg.config.software_params.dup_concordance_cutoff,
    output:
        known="sample_concordance/KnownReplicates.csv",
        known_qc="sample_concordance/InternalQcKnown.csv",
        known_study="sample_concordance/StudySampleKnown.csv",
        unknown="sample_concordance/UnknownReplicates.csv",
    script:
        "../scripts/known_concordant_samples.py"
