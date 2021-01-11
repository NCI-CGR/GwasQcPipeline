"""Reporting modules

This modules contains all things reporting.
"""
import pandas as pd


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
        all_samples="sample_qc_report/all_samples.csv",
    script:
        "../scripts/sample_qc_report.py"


rule lims_upload:
    input:
        rules.sample_qc_report.output.all_samples,
    output:
        "sample_qc_report/lims_upload.csv",
    run:
        (
            pd.read_csv(input[0])
            .rename({"Call_Rate_Initial": "Call Rate"}, axis=1)
            .reindex(
                [
                    "SR_Subject_ID",
                    "LIMS_Individual_ID",
                    "Sample_ID",
                    "Project-Sample ID",
                    "Call Rate",
                    "Low Call Rate",
                    "Contaminated",
                    "Sex Discordant",
                    "Expected Replicate Discordance",
                    "Unexpected Replicate",
                ],
                axis=1,
            )
            .to_csv(output[0], index=False)
        )


rule identifiler_needed:
    input:
        rules.sample_qc_report.output.all_samples,
    output:
        "sample_qc_report/identifiler_needed.csv",
    run:
        (
            pd.read_csv(input[0])
            .query("Identifiler_Needed")
            .reindex(
                [
                    "Sample_ID",
                    "LIMSSample_ID",
                    "Project",
                    "Project-Sample ID",
                    "SR_Subject_ID",
                    "LIMS_Individual_ID",
                    "Identifiler_Reason",
                ],
                axis=1,
            )
            .to_csv(output[0], index=False)
        )


rule sample_qc_report_summary_stats:
    input:
        rules.sample_qc_report.output.all_samples,
    output:
        "sample_qc_report/summary_stats.txt",
    script:
        "../scripts/sample_qc_report_summary_stats.py"
