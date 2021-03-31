import pandas as pd

from cgr_gwas_qc.reporting import REPORT_NAME_MAPPER


################################################################################
# Files For Lab
################################################################################
rule lab_sample_level_qc_report:
    input:
        "sample_level/sample_qc.csv",
    output:
        "files_for_lab/{deliver_prefix}all_sample_qc{deliver_suffix}.csv",
    shell:
        "cp {input[0]} {output[0]}"


rule lab_lims_upload:
    input:
        "sample_level/sample_qc.csv",
    output:
        "files_for_lab/{deliver_prefix}LimsUpload{deliver_suffix}.csv",
    run:
        (
            pd.read_csv(input[0])
            .rename(REPORT_NAME_MAPPER, axis=1)
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
        "sample_level/sample_qc.csv",
    output:
        "files_for_lab/{deliver_prefix}Identifiler{deliver_suffix}.csv",
    run:
        (
            pd.read_csv(input[0])
            .query("identifiler_needed")
            .rename(REPORT_NAME_MAPPER, axis=1)
            .reindex(
                [
                    "Sample_ID",
                    "LIMSSample_ID",
                    "Project",
                    "Project-Sample ID",
                    "SR_Subject_ID",
                    "LIMS_Individual_ID",
                    "identifiler_reason",
                ],
                axis=1,
            )
            .to_csv(output[0], index=False)
        )


rule lab_known_replicates:
    input:
        "sample_level/concordance/KnownReplicates.csv",
    output:
        "files_for_lab/{deliver_prefix}KnownReplicates{deliver_suffix}.csv",
    shell:
        "cp {input[0]} {output[0]}"


rule lab_unknown_replicates:
    input:
        "sample_level/concordance/UnknownReplicates.csv",
    output:
        "files_for_lab/{deliver_prefix}UnknownReplicates{deliver_suffix}.csv",
    shell:
        "cp {input[0]} {output[0]}"


################################################################################
# Deliverables
################################################################################
rule deliver_manifest:
    input:
        cfg.sample_sheet_file.as_posix(),
    output:
        "deliver/{deliver_prefix}AnalysisManifest{deliver_suffix}.csv",
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
        "sample_level/sample_qc.csv",
    output:
        "deliver/SampleUsedforSubject.csv",
    run:
        qc = pd.read_csv(input[0]).query("not is_internal_control")  # exclude internal controls

        (
            qc.query("is_subject_representative")
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


################################################################################
# QC Report
################################################################################
rule qc_report:
    input:
        sample_sheet_csv=cfg.sample_sheet_file,
        snp_qc_csv="sample_level/snp_qc.csv",
        sample_qc_csv="sample_level/sample_qc.csv",
        population_qc_csv="population_level/population_qc.csv",
        control_replicates_csv="sample_level/concordance/InternalQcKnown.csv",
        study_replicates_csv="sample_level/concordance/StudySampleKnown.csv",
        unexpected_replicates_csv="sample_level/concordance/UnknownReplicates.csv",
        call_rate_png="sample_level/call_rate.png",
        chrx_inbreeding_png="sample_level/chrx_inbreeding.png",
        ancestry_png="sample_level/ancestry.png",
        autosomal_heterozygosity_png_dir="population_level/autosomal_heterozygosity_plots",
        pca_png_dir="population_level/pca_plots",
        hwe_png_dir="population_level/hwe_plots",
    params:
        config=cfg.config,
    output:
        "deliver/qc_report.md",
    script:
        "../scripts/qc_report.py"


rule export_qc_report:
    input:
        rules.qc_report.output[0],
    output:
        "deliver/{deliver_prefix}QC_Report{deliver_suffix}.{ext}",
    wildcard_constraints:
        ext="docx|pdf|html",
    conda:
        cfg.conda("pandoc.yml")
    shell:
        "pandoc --toc -s {input} -o {output[0]}"
