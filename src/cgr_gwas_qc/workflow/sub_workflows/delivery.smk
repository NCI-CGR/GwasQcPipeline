from cgr_gwas_qc import load_config
cfg = load_config()
output_pattern = cfg.config.user_files.output_pattern
manifest_file = cfg.config.sample_sheet

localrules:
    all_delivery,
    lab_sample_level_qc_report,
    lab_lims_upload,
    lab_identifiler_needed,
    lab_known_replicates,
    lab_unknown_replicates,
    deliver_manifest,
    deliver_hwp,
    deliver_original_sample_data,
    deliver_subject_data,
    deliver_subject_list,
    qc_report,
    qc_report_docx,
    qc_report_xlsx,


wildcard_constraints:
    deliver_prefix=".*",
    deliver_suffix=".*",


################################################################################
# Delivery Targets
################################################################################
targets = [
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
]

if cfg.config.workflow_params.lims_upload:
    # The CGEMs/CCAD cluster has a cron job running that looks for this file in
    # the root run directory. If it is there then it will automatically upload
    # to the LIMs system. This is only useful on CGEMs/CCAD.
    targets.append(output_pattern.format(prefix=".", file_type="LimsUpload", ext="csv"))


rule all_delivery:
    input:
        targets,


################################################################################
# Workflow Rules
################################################################################
# -------------------------------------------------------------------------------
# Files For Lab
# -------------------------------------------------------------------------------
rule lab_sample_level_qc_report:
    input:
        "sample_level/sample_qc.csv",
    output:
        "files_for_lab/{deliver_prefix}all_sample_qc{deliver_suffix}.csv",
    shell:
        "cp {input[0]} {output[0]}"


rule lab_lims_upload:
    input:
        sample_sheet_csv="cgr_sample_sheet.csv",
        sample_qc_csv="sample_level/sample_qc.csv",
    output:
        "files_for_lab/{deliver_prefix}LimsUpload{deliver_suffix}.csv",
    script:
        "../scripts/lab_lims_upload.py"


if cfg.config.workflow_params.lims_upload:

    ruleorder: lab_lims_upload > lims_upload

    rule lims_upload:
        input:
            rules.lab_lims_upload.output[0],
        output:
            "{deliver_prefix}LimsUpload{deliver_suffix}.csv",
        shell:
            "cp {input[0]} {output[0]}"


rule lab_identifiler_needed:
    input:
        sample_sheet_csv="cgr_sample_sheet.csv",
        sample_qc_csv="sample_level/sample_qc.csv",
    output:
        "files_for_lab/{deliver_prefix}Identifiler{deliver_suffix}.csv",
    script:
        "../scripts/lab_identifiler_needed.py"


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


# -------------------------------------------------------------------------------
# Deliverables
# -------------------------------------------------------------------------------
rule deliver_manifest:
    input:
        cfg.sample_sheet_file.as_posix(),
        manifest = manifest_file,
    output:
        "delivery/{deliver_prefix}AnalysisManifest{deliver_suffix}.csv",
    shell:
        """
        cp {input[0]} {output[0]}
        file=$(echo {input[1]} | rev | cut -d"/" -f -1 | rev)
        cp {input[1]} delivery/$file
        """
        


rule deliver_hwp:
    input:
        "subject_level/.control_plots.done",
    output:
        "delivery/HWP.zip",
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
        bed="sample_level/samples.bed",
        bim="sample_level/samples.bim",
        fam="sample_level/samples.fam",
    output:
        bed="delivery/samples.bed",
        bim="delivery/samples.bim",
        fam="delivery/samples.fam",
    shell:
        "cp {input.bed} {output.bed} && cp {input.bim} {output.bim} && cp {input.fam} {output.fam}"


rule deliver_subject_data:
    input:
        bed="subject_level/samples.bed",
        bim="subject_level/samples.bim",
        fam="subject_level/samples.fam",
    output:
        bed="delivery/subjects.bed",
        bim="delivery/subjects.bim",
        fam="delivery/subjects.fam",
    shell:
        "cp {input.bed} {output.bed} && cp {input.bim} {output.bim} && cp {input.fam} {output.fam}"


rule deliver_subject_list:
    input:
        "sample_level/sample_qc.csv",
    output:
        "delivery/SampleUsedforSubject.csv",
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
        snp_qc_csv="sample_level/snp_qc.csv",
        sample_qc_csv="sample_level/sample_qc.csv",
        subject_qc_csv="subject_level/subject_qc.csv",
        population_qc_csv="subject_level/population_qc.csv",
        control_replicates_csv="sample_level/concordance/InternalQcKnown.csv",
        study_replicates_csv="sample_level/concordance/StudySampleKnown.csv",
        unexpected_replicates_csv="sample_level/concordance/UnknownReplicates.csv",
        call_rate_png="sample_level/call_rate.png",
        chrx_inbreeding_png="sample_level/chrx_inbreeding.png",
        ancestry_png="sample_level/ancestry.png",
        _population_plots="subject_level/.population_plots.done",
        _control_plots="subject_level/.control_plots.done",
    params:
        config=cfg.config,
        autosomal_heterozygosity_png_dir="subject_level/autosomal_heterozygosity_plots",
        pca_png_dir="subject_level/pca_plots",
        hwe_png_dir="subject_level/hwe_plots",
    output:
        "delivery/qc_report.md",
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
    shell:
        """
        pandoc --reference-doc {params.template} --toc -s {input} -o {output[0]}
        #ml  plink/1.9
	#plink \
	#--bfile sample_level/samples \
	#--pheno subject_level/gwas.txt \
	#--ci 0.95 \
	#--assoc \
	#--out delivery/gwas_subjects \
	#--memory 20000 \
	#--threads 2 \
	#--allow-no-sex
	#ml -plink/1.9

        """


rule qc_report_xlsx:
    input:
        sample_sheet_csv="cgr_sample_sheet.csv",
        sample_concordance_csv="sample_level/concordance/summary.csv",
        sample_qc_csv="sample_level/sample_qc.csv",
        subject_qc_csv="subject_level/subject_qc.csv",
        population_concordance_csv="subject_level/concordance.csv",
        population_qc_csv="subject_level/population_qc.csv",
        graf="sample_level/ancestry/graf_populations.txt",
    output:
        "delivery/{deliver_prefix}QC_Report{deliver_suffix}.xlsx",
    script:
        "../scripts/qc_report_table.py"
