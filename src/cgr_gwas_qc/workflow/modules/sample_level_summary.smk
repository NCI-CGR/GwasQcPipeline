################################################################################
# Remove Contaminated
################################################################################
rule Sample_IDs_above_contam_threshold:
    """Creates a list of contaminated Sample_IDs.

    Compares the %Mix from ``verifyIDintensity`` to the contamination
    threshold from the config. If %Mix is greater than this threshold
    Sample_IDs are saved to a list.
    """
    input:
        rules.agg_contamination_test.output[0],
    params:
        contam_threshold=cfg.config.software_params.contam_threshold,
    output:
        "sample_contamination/contaminated_samples.txt",
    run:
        sample_ids = (
            pd.read_csv(input[0]).query(f"`%Mix` > {params.contam_threshold}").Sample_ID.values
        )
        with open(output[0], "w") as fh:
            for sample_id in sample_ids:
                fh.write("{} {}\n".format(sample_id, sample_id))


rule remove_contaminated_samples:
    input:
        bed="plink_filter_call_rate_2/samples.bed",
        bim="plink_filter_call_rate_2/samples.bim",
        fam="plink_filter_call_rate_2/samples.fam",
        contaminated_Sample_IDs=rules.Sample_IDs_above_contam_threshold.output[0],
    params:
        in_prefix="plink_filter_call_rate_2/samples",
        out_prefix="sample_contamination/not_contaminated/samples",
    output:
        bed="sample_contamination/not_contaminated/samples.bed",
        bim="sample_contamination/not_contaminated/samples.bim",
        fam="sample_contamination/not_contaminated/samples.fam",
        nosex="sample_contamination/not_contaminated/samples.nosex",
    log:
        "sample_contamination/not_contaminated/samples.log",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    threads: 2
    resources:
        mem=10000,
    shell:
        "plink "
        "--bfile {params.in_prefix} "
        "--remove {input.contaminated_Sample_IDs} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem} "
        "--out {params.out_prefix}"
