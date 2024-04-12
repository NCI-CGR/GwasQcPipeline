import pandas as pd
import seaborn as sns

CASE_CONTROL_DTYPE = pd.CategoricalDtype(categories=["Case", "Control", "QC", "Unknown"])
CASE_CONTROL_COLORS = ["#f7022a", "#3e82fc", "gray", "gold"]  # red  # blue  # gray  #gold

# Assign labels to colors for plotting consistency
CASE_CONTROL_LABEL_COLORS = {'Case': CASE_CONTROL_COLORS[0], 'Control': CASE_CONTROL_COLORS[1], 'QC': CASE_CONTROL_COLORS[2], 'Unknown': CASE_CONTROL_COLORS[3]}

SEX_DTYPE = pd.CategoricalDtype(categories=["F", "M", "U"])

# Mapping current column names to names from the legacy workflow to maintain
# consistency in deliverables.
REPORT_NAME_MAPPER = {
    "case_control": "Case/Control_Status",
    "concordance": "Concordance",
    "expected_sex": "Expected_Sex",
    "identifiler_needed": "Identifiler_Needed",
    "identifiler_reason": "Identifiler Reason",
    "is_call_rate_filtered": "Low Call Rate",
    "is_contaminated": "Contaminated",
    "is_cr1_filtered": "Call_Rate_1_filter",
    "is_cr2_filtered": "Call_Rate_2_filter",
    "is_internal_control": "Internal Control",
    "is_expected_replicate": "Expected Replicate",
    "is_discordant_replicate": "Expected Replicate Discordance",
    "is_sample_exclusion": "Sample Excluded from QC",
    "is_sex_discordant": "Sex Discordant",
    "is_subject_representative": "Sample Used in Subject Analysis",
    "is_unexpected_replicate": "Unexpected Replicate",
    "num_samples_per_subject": "Number Samples Per Subject",
    "predicted_sex": "Predicted_Sex",
    "QC_Family_ID": "Fam_ID",
    "replicate_ids": "Replicate IDs",
    "subject_dropped_from_study": "Subject Removed",
    "X_inbreeding_coefficient": "ChrX_Inbreed_estimate",
}
