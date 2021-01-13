
rule sample_concordance:
    """Summarize sample concordance using IBS/IBD.

    Calculates the proportion of shared homozygous markers (IBS2 / (IBS0 + IBS1 + IBS2)) as a
    measure of sample concordance. Then outputs concordance measures for samples that are known
    to be replicates and samples that are thought to be unrelated/independent with a concordance
    > dup_concordance_cutoff (currently 0.95).
    """
    input:
        sample_sheet=cfg.sample_sheet_file,
        imiss="sample_level/call_rate_2/samples.imiss",
        ibd="sample_level/call_rate_2/samples_maf{maf}_ld{ld}.genome".format(
            maf=cfg.config.software_params.maf_for_ibd, ld=cfg.config.software_params.ld_prune_r2
        ),
    params:
        subject_id_override=cfg.config.workflow_params.subject_id_to_use,
        concordance_threshold=cfg.config.software_params.dup_concordance_cutoff,
    output:
        known="sample_level/concordance/KnownReplicates.csv",
        known_qc="sample_level/concordance/InternalQcKnown.csv",
        known_study="sample_level/concordance/StudySampleKnown.csv",
        unknown="sample_level/concordance/UnknownReplicates.csv",
    script:
        "../scripts/known_concordant_samples.py"
