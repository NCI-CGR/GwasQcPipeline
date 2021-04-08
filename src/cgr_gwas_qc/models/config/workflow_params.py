import warnings
from typing import Optional

from pydantic import BaseModel, Field, validator


class WorkflowParams(BaseModel):
    """This set of parameters control what parts and how the workflow is run."""

    subject_id_column: Optional[str] = Field(
        None,
        description="[Deprecated] The name of the column in the sample sheet which identifies unique subjects.",
    )
    expected_sex_column: str = Field(
        "Expected_Sex", description="Column in the sample sheet that describes the expected sex."
    )

    remove_contam: bool = Field(
        True,
        description="Contaminated samples will be removed prior to sample to subject transformation.",
    )
    remove_sex_discordant: bool = Field(
        True,
        description="Sex discordant samples will be removed prior to sample to subject transformation.",
    )
    remove_rep_discordant: bool = Field(
        True,
        description="Discordant replicates will be pruned prior to sample to subject transformation.",
    )
    remove_unexpected_rep: bool = Field(
        True,
        description="Unexpected replicates will be pruned prior to sample to subject transformation.",
    )

    minimum_pop_subjects: int = Field(
        50, description="Minimum number of subjects required in order to analyze a population", gt=0
    )
    control_hwp_threshold: int = Field(
        50,
        description="Minimum number of controls (in a population) required for HWE estimation",
        gt=0,
    )

    @validator("subject_id_column")
    def validate_subject_id_column(cls, v):
        if v is None:
            return v

        warnings.warn(
            "subject_id_column is deprecated, add this to the Group_By column in the sample sheet.",
            DeprecationWarning,
        )

        return v
