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
        lims=temp("qc_table_for_LimsUpload.csv"),
    script:
        "../scripts/sample_qc_report.py"


def _rename_lims_upload():
    ss_file = cfg.sample_sheet_file.stem
    if "_AnalysisManifest_" in ss_file:
        out_name, sheet_date = ss_file.split("_AnalysisManifest_")
        return f"{out_name}_LimsUpload_{sheet_date}.csv"
    return "all_sample_qc_LimsUpload.csv"


rule rename_lims_upload:
    input:
        rules.sample_qc_report.output.lims,
    output:
        _rename_lims_upload(),
    shell:
        "cp {input[0]} {output[0]}"


rule summary_stats:
    input:
        "all_sample_qc.csv",
    output:
        "summary_stats.txt",
    script:
        "../scripts/sample_qc_report_summary_stats.py"
