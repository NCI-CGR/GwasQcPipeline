import pandas as pd


def use_lims_name(pth, file_type, ext="csv"):
    if "AnalysisManifest" in cfg.sample_sheet_file.stem:
        lims_name = cfg.sample_sheet_file.stem.replace("AnalysisManifest", file_type)
        return f"{pth}/{lims_name}.{ext}"
    return f"{pth}/{file_type}.{ext}"


################################################################################
# Files For Lab
################################################################################
rule lab_sample_level_qc_report:
    input:
        "sample_level/qc_summary.csv",
    output:
        use_lims_name("files_for_lab", "all_sample_qc"),
    shell:
        "cp {input[0]} {output[0]}"


rule lab_lims_upload:
    input:
        "sample_level/qc_summary.csv",
    output:
        use_lims_name("files_for_lab", "LimsUpload"),
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


rule lab_identifiler_needed:
    input:
        "sample_level/qc_summary.csv",
    output:
        use_lims_name("files_for_lab", "Identifiler"),
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


rule lab_known_replicates:
    input:
        "sample_level/concordance/KnownReplicates.csv",
    output:
        use_lims_name("files_for_lab", "KnownReplicates"),
    shell:
        "cp {input[0]} {output[0]}"


rule lab_unknown_replicates:
    input:
        "sample_level/concordance/UnknownReplicates.csv",
    output:
        use_lims_name("files_for_lab", "UnknownReplicates"),
    shell:
        "cp {input[0]} {output[0]}"


################################################################################
# Deliverables
################################################################################
rule deliver_manifest:
    input:
        cfg.sample_sheet_file.as_posix(),
    output:
        use_lims_name("deliver", "AnalysisManifest"),
    shell:
        "cp {input[0]} {output[0]}"


rule deliver_hwp:
    input:
        "population_level/controls.done",
    output:
        "deliver/HWP.zip",
    shell:
        """
        ODIR=$(dirname {output[0]})
        mkdir -p $ODIR/HWP \
        && for FILE in $(cat {input[0]}); do cp $FILE $ODIR/HWP/; done \
        && cd $ODIR \
        && zip -r $(basename {output[0]}) ./HWP \
        && rm -r HWP
        """


rule deliver_original_sample_data:
    input:
        bed="sample_level/samples.bed",
        bim="sample_level/samples.bim",
        fam="sample_level/samples.fam",
    output:
        bed="deliver/samples.bed",
        bim="deliver/samples.bim",
        fam="deliver/samples.fam",
    shell:
        "cp {input.bed} {output.bed} && cp {input.bim} {output.bim} && cp {input.fam} {output.fam}"


rule deliver_subject_data:
    input:
        bed="subject_level/samples.bed",
        bim="subject_level/samples.bim",
        fam="subject_level/samples.fam",
    output:
        bed="deliver/subjects.bed",
        bim="deliver/subjects.bim",
        fam="deliver/subjects.fam",
    shell:
        "cp {input.bed} {output.bed} && cp {input.bim} {output.bim} && cp {input.fam} {output.fam}"


rule deliver_subject_list:
    input:
        "sample_level/qc_summary.csv",
    output:
        "deliver/SampleUsedforSubject.csv",
    run:
        qc = pd.read_csv(input[0]).query("not Internal_Control")  # exclude internal controls

        (
            qc.query("Subject_Representative")
            .reindex(["Group_By_Subject_ID", "Sample_ID"], axis=1)
            # For groups without representative set Sample_ID to NA
            .set_index("Group_By_Subject_ID")
            .reindex(qc.Group_By_Subject_ID.unique())
            .fillna("NA")
            .reset_index()
            # Clean-up
            .rename({"Group_By_Subject_ID": "Subject_ID"}, axis=1)
            .to_csv(output[0], index=False)
        )


rule deliver_readme:
    output:
        "deliver/README",


rule deliver_summary_writeup:
    output:
        use_lims_name("deliver", "QC_Report", "docx"),


rule deliver_summary_table:
    output:
        use_lims_name("deliver", "QC_Report", "xlsx"),
