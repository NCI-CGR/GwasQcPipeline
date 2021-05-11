from cgr_gwas_qc import load_config

cfg = load_config()
output_pattern = cfg.config.user_files.output_pattern


wildcard_constraints:
    deliver_prefix=".*",
    deliver_suffix=".*",


localrules:
    all_delivery,


################################################################################
# All Targets
################################################################################
rule all_delivery:
    input:
        "delivery/samples.bed",
        "delivery/samples.bim",
        "delivery/samples.fam",
        "delivery/SampleUsedforSubject.csv",
        "delivery/subjects.bed",
        "delivery/subjects.bim",
        "delivery/subjects.fam",
        "delivery/HWP.zip",
        output_pattern.format(prefix="files_for_lab", file_type="all_sample_qc", ext="csv"),
        output_pattern.format(prefix="files_for_lab", file_type="LimsUpload", ext="csv"),
        output_pattern.format(prefix="files_for_lab", file_type="Identifiler", ext="csv"),
        output_pattern.format(prefix="files_for_lab", file_type="KnownReplicates", ext="csv"),
        output_pattern.format(prefix="files_for_lab", file_type="UnknownReplicates", ext="csv"),
        output_pattern.format(prefix="delivery", file_type="AnalysisManifest", ext="csv"),
        output_pattern.format(prefix="delivery", file_type="QC_Report", ext="docx"),
        output_pattern.format(prefix="delivery", file_type="QC_Report", ext="xlsx"),


################################################################################
# Imports
################################################################################
subworkflow entry_points:
    snakefile:
        cfg.subworkflow("entry_points")
    workdir:
        cfg.root.as_posix()


subworkflow sample_qc:
    snakefile:
        cfg.subworkflow("sample_qc")
    workdir:
        cfg.root.as_posix()


subworkflow subject_qc:
    snakefile:
        cfg.subworkflow("subject_qc")
    workdir:
        cfg.root.as_posix()


################################################################################
# Workflow Rules
################################################################################
# -------------------------------------------------------------------------------
# Files For Lab
# -------------------------------------------------------------------------------
rule lab_sample_level_qc_report:
    input:
        sample_qc("sample_level/sample_qc.csv"),
    output:
        "files_for_lab/{deliver_prefix}all_sample_qc{deliver_suffix}.csv",
    group:
        "delivery"
    shell:
        "cp {input[0]} {output[0]}"


rule lab_lims_upload:
    input:
        sample_sheet_csv="cgr_sample_sheet.csv",
        sample_qc_csv=sample_qc("sample_level/sample_qc.csv"),
    output:
        "files_for_lab/{deliver_prefix}LimsUpload{deliver_suffix}.csv",
    group:
        "delivery"
    script:
        "../scripts/lab_lims_upload.py"


rule lab_identifiler_needed:
    input:
        sample_sheet_csv="cgr_sample_sheet.csv",
        sample_qc_csv=sample_qc("sample_level/sample_qc.csv"),
    output:
        "files_for_lab/{deliver_prefix}Identifiler{deliver_suffix}.csv",
    group:
        "delivery"
    script:
        "../scripts/lab_identifiler_needed.py"


rule lab_known_replicates:
    input:
        sample_qc("sample_level/concordance/KnownReplicates.csv"),
    output:
        "files_for_lab/{deliver_prefix}KnownReplicates{deliver_suffix}.csv",
    group:
        "delivery"
    shell:
        "cp {input[0]} {output[0]}"


rule lab_unknown_replicates:
    input:
        sample_qc("sample_level/concordance/UnknownReplicates.csv"),
    output:
        "files_for_lab/{deliver_prefix}UnknownReplicates{deliver_suffix}.csv",
    group:
        "delivery"
    shell:
        "cp {input[0]} {output[0]}"


# -------------------------------------------------------------------------------
# Deliverables
# -------------------------------------------------------------------------------
rule deliver_manifest:
    input:
        cfg.sample_sheet_file.as_posix(),
    output:
        "delivery/{deliver_prefix}AnalysisManifest{deliver_suffix}.csv",
    group:
        "delivery"
    shell:
        "cp {input[0]} {output[0]}"


rule deliver_hwp:
    input:
        subject_qc("subject_level/.control_plots.done"),
    output:
        "delivery/HWP.zip",
    group:
        "delivery"
    shell:
        """
        ODIR=$(dirname {output[0]})
        mkdir -p $ODIR/HWP
        for FILE in $(find ./subject_level -type f -name '*.hwe'); do
            cp $FILE $ODIR/HWP/
        done
        cd $ODIR
        zip -r $(basename {output[0]}) ./HWP
        rm -r HWP
        """


rule deliver_original_sample_data:
    input:
        bed=entry_points("sample_level/samples.bed"),
        bim=entry_points("sample_level/samples.bim"),
        fam=entry_points("sample_level/samples.fam"),
    output:
        bed="delivery/samples.bed",
        bim="delivery/samples.bim",
        fam="delivery/samples.fam",
    group:
        "delivery"
    shell:
        "cp {input.bed} {output.bed} && cp {input.bim} {output.bim} && cp {input.fam} {output.fam}"


rule deliver_subject_data:
    input:
        bed=subject_qc("subject_level/samples.bed"),
        bim=subject_qc("subject_level/samples.bim"),
        fam=subject_qc("subject_level/samples.fam"),
    output:
        bed="delivery/subjects.bed",
        bim="delivery/subjects.bim",
        fam="delivery/subjects.fam",
    group:
        "delivery"
    shell:
        "cp {input.bed} {output.bed} && cp {input.bim} {output.bim} && cp {input.fam} {output.fam}"


rule deliver_subject_list:
    input:
        sample_qc("sample_level/sample_qc.csv"),
    output:
        "delivery/SampleUsedforSubject.csv",
    group:
        "delivery"
    run:
        import pandas as pd

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


# -------------------------------------------------------------------------------
# QC Report
# -------------------------------------------------------------------------------
rule qc_report:
    input:
        sample_sheet_csv="cgr_sample_sheet.csv",
        snp_qc_csv=sample_qc("sample_level/snp_qc.csv"),
        sample_qc_csv=sample_qc("sample_level/sample_qc.csv"),
        subject_qc_csv=subject_qc("subject_level/subject_qc.csv"),
        population_qc_csv=subject_qc("subject_level/population_qc.csv"),
        control_replicates_csv=sample_qc("sample_level/concordance/InternalQcKnown.csv"),
        study_replicates_csv=sample_qc("sample_level/concordance/StudySampleKnown.csv"),
        unexpected_replicates_csv=sample_qc("sample_level/concordance/UnknownReplicates.csv"),
        call_rate_png=sample_qc("sample_level/call_rate.png"),
        chrx_inbreeding_png=sample_qc("sample_level/chrx_inbreeding.png"),
        ancestry_png=sample_qc("sample_level/ancestry.png"),
        _population_plots=subject_qc("subject_level/.population_plots.done"),
        _control_plots=subject_qc("subject_level/.population_plots.done"),
    params:
        config=cfg.config,
        autosomal_heterozygosity_png_dir="subject_level/autosomal_heterozygosity_plots",
        pca_png_dir="subject_level/pca_plots",
        hwe_png_dir="subject_level/hwe_plots",
    output:
        "delivery/qc_report.md",
    group:
        "delivery"
    script:
        "../scripts/qc_report.py"


rule qc_report_docx:
    input:
        rules.qc_report.output[0],
    params:
        template=cfg.docx_template,
    output:
        "delivery/{deliver_prefix}QC_Report{deliver_suffix}.docx",
    conda:
        cfg.conda("pandoc")
    group:
        "delivery"
    shell:
        "pandoc --reference-doc {params.template} --toc -s {input} -o {output[0]}"


rule qc_report_xlsx:
    input:
        sample_sheet_csv="cgr_sample_sheet.csv",
        sample_concordance_csv=sample_qc("sample_level/concordance/summary.csv"),
        sample_qc_csv=sample_qc("sample_level/sample_qc.csv"),
        subject_qc_csv=subject_qc("subject_level/subject_qc.csv"),
        population_concordance_csv=subject_qc("subject_level/concordance.csv"),
        population_qc_csv=subject_qc("subject_level/population_qc.csv"),
        graf=sample_qc("sample_level/ancestry/graf_populations.txt"),
    output:
        "delivery/{deliver_prefix}QC_Report{deliver_suffix}.xlsx",
    group:
        "delivery"
    script:
        "../scripts/qc_report_table.py"
