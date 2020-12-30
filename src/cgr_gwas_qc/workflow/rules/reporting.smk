"""Reporting modules

This modules contains all things reporting.
"""


def sample_qc_report_inputs(wildcards):
    inputs = {
        "imiss_start": "plink_start/samples.imiss",
        "imiss_cr1": "plink_filter_call_rate_1/samples.imiss",
        "imiss_cr2": "plink_filter_call_rate_2/samples.imiss",
        "sexcheck_cr1": "plink_filter_call_rate_1/samples.sexcheck",
        "ancestry": "ancestry/graf_ancestry_calls.txt",
        "known_replicates": "sample_filters/concordance/KnownReplicates.csv",
        "unknown_replicates": "sample_filters/concordance/UnknownReplicates.csv",
    }

    if (
        cfg.config.user_files.idat_pattern
        and cfg.config.user_files.gtc_pattern
        and cfg.config.workflow_params.remove_contam
    ):
        inputs["contam"] = "sample_filters/agg_contamination_test.csv"
        inputs["intensity"] = "sample_filters/agg_median_idat_intensity.csv"

    return inputs


rule sample_qc_report:
    input:
        unpack(sample_qc_report_inputs),
    output:
        all_qc="all_sample_qc.csv",
        lims="{out_name}_LimsUpload_{sample_sheet_date}.csv",
    script:
        "../scripts/sample_qc_report.py"
