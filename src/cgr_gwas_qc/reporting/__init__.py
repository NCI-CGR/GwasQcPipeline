import pandas as pd
from jinja2 import Environment, PackageLoader

env = Environment(
    loader=PackageLoader("cgr_gwas_qc", "reporting"),
    trim_blocks=True,
    lstrip_blocks=True,
    keep_trailing_newline=True,
)

CASE_CONTROL_DTYPE = pd.CategoricalDtype(categories=["Case", "Control", "QC"])
SEX_DTYPE = pd.CategoricalDtype(categories=["M", "F", "U"])

# Mapping current column names to names from the legacy workflow to maintain
# consistency in deliverables.
REPORT_NAME_MAPPER = {
    "expected_sex": "Expected_Sex",
    "predicted_sex": "Predicted_Sex",
    "sex_discordant": "Sex Discordant",
    "X_inbreeding_coefficient": "ChrX_Inbreed_estimate",
    "is_call_rate_filtered": "Low Call Rate",
    "idats_exist": "IdatsInProjectDir",
}
