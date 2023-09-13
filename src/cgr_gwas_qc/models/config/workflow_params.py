from pydantic import BaseModel, Field


class WorkflowParams(BaseModel):
    """This set of parameters control what parts and how the workflow is run.

    .. code-block:: yaml

        workflow_params:
            subject_id_column: Group_By
            expected_sex_column: Expected_Sex
            sex_chr_included: true
            case_control_column: Case/Control_Status
            remove_contam: true
            remove_rep_discordant: true
            minimum_pop_subjects: 50
            control_hwp_threshold: 50
            lims_upload: false
            case_control_gwas: false
    """

    subject_id_column: str = Field(
        "Group_By",
        description="The column in your sample sheet that contains the subject ID. "
        "We expect that there may be multiple samples (rows) that have the sample subject ID. "
        "If there are multiple columns that you want use for subject ID, you can create a special column named ``Group_By``. "
        "This would be a column that contains the column name to use as subject ID for that given samples (row).",
    )

    expected_sex_column: str = Field(
        "Expected_Sex",
        description="The name of the column in the sample sheet which identifies expected sex of samples. "
        "Allowed values in this columns are [``F`` | ``M`` | ``U``]",
    )

    sex_chr_included: bool = Field(
        True,
        description="True if the sex chromosome is included in the microarray and a sex concordance check can be performed.",
    )

    case_control_column: str = Field(
        "Case/Control_Status",
        description="The name of the colun in the sample sheet which identifies Case/Control status. "
        "Allowed values in this column are [``Control`` | ``Case`` | ``QC`` | ``Unknown``].",
    )

    remove_contam: bool = Field(
        True,
        description="True if you want to remove contaminated samples before running the subject level QC.",
    )

    remove_rep_discordant: bool = Field(
        True,
        description="True if you want to remove discordant replicates before running the subject level QC.",
    )

    minimum_pop_subjects: int = Field(
        50,
        gt=0,
        description="The minimum number of samples needed to use a population for population level QC (PCA, Autosomal Heterozygosity).",
    )

    control_hwp_threshold: int = Field(
        50,
        gt=0,
        description="The minimum number of control samples needed to use a population for population level QC (HWE).",
    )

    lims_upload: bool = Field(
        False,
        description="For ``CGEMS/CCAD`` use only, will place a copy of the LimsUpload file in the root directory.",
    )

    case_control_gwas: bool = Field(
        False,
        description="A plink logistic regression gwas will be performed with case_control phenotype.",
    )

    @staticmethod
    def schema_rst():
        """Tweak schema for rendering in Sphinx."""
        import copy
        import json

        content = copy.deepcopy(WorkflowParams.schema())
        content["title"] = "Workflow Parameters"

        return json.dumps(content, indent=2)
