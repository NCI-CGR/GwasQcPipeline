from pydantic import BaseModel, Field


class WorkflowParams(BaseModel):
    """This set of parameters control what parts and how the workflow is run."""

    subject_id_column: str = Field(
        "Group_By",
        description="The name of the column in the sample sheet which identifies unique subjects.",
    )

    expected_sex_column: str = Field(
        "Expected_Sex",
        description="The name of the column in the sample sheet which identifies expected sex of samples.",
    )

    case_control_column: str = Field(
        "Case/Control_Status",
        description="The name of the column in the sample sheet which identifies Case/Control status.",
    )

    remove_contam: bool = Field(
        True,
        description="Contaminated samples will be removed prior to sample to subject transformation.",
    )
    remove_rep_discordant: bool = Field(
        True,
        description="Discordant replicates will be pruned prior to sample to subject transformation.",
    )

    minimum_pop_subjects: int = Field(
        50, description="Minimum number of subjects required in order to analyze a population", gt=0
    )
    control_hwp_threshold: int = Field(
        50,
        description="Minimum number of controls (in a population) required for HWE estimation",
        gt=0,
    )

    lims_upload: bool = Field(
        False,
        description="Create summary file for automatic upload to the CGEMs/CCAD LIMs system.",
    )
