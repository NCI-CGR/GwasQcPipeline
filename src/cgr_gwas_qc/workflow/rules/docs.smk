import pandas as pd


rule subset_all_qc:
    input:
        "sample_qc_report/all_samples.csv",
    output:
        "word_doc/all_sample_qc.csv",
    run:
        (
            pd.read_csv(input[0])
            .drop(
                [
                    "Count_of_QC_Issue",
                    "Count_of_SR_SubjectID",
                    "Current_Subject_Status",
                    "Expected Replicate Discordance",
                    "Identifiler_Needed",
                    "SR",
                    "Subject_Notes",
                ],
                axis=1,
                errors="ignore",
            )
            .to_csv(output[0], index=False)
        )
