import pandas as pd

CASE_CONTROL_DTYPE = pd.CategoricalDtype(categories=["Case", "Control", "QC", "Unknown"])
CASE_CONTROL_COLORS = ["#f7022a", "#3e82fc", "gray", "#1bfc06"]  # red  # blue  # gray  # green

SEX_DTYPE = pd.CategoricalDtype(categories=["F", "M", "U"])

# Mapping current column names to names from the legacy workflow to maintain
# consistency in deliverables.
REPORT_NAME_MAPPER = {
    "X_inbreeding_coefficient": "ChrX_Inbreed_estimate",
    "expected_sex": "Expected_Sex",
    "predicted_sex": "Predicted_Sex",
    "idats_exist": "IdatsInProjectDir",
    "identifiler_needed": "Identifiler_Needed",
    "is_call_rate_filtered": "Low Call Rate",
    "is_contaminated": "Contaminated",
    "is_replicate_discordant": "Expected Replicate Discordance",
    "is_sex_discordant": "Sex Discordant",
    "is_unexpected_replicate": "Unexpected Replicate",
}
