import time
from typing import Optional

from pydantic import BaseModel, Field

# setting start time to define unique lim file name
timestr = time.strftime("%Y%m%d%H%M%S")


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
            lims_upload: true
            lims_output_dir: /DCEG/CGF/Laboratory/LIMS/drop-box-prod/gwas_primaryqc/
            case_control_gwas: false
            max_time_hr:
            max_mem_mb:
            time_start:
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
        True,
        description="For ``CGEMS/CCAD`` use only, will place a copy of the LimsUpload file in the root directory.",
    )

    lims_output_dir: str = Field(
        "/DCEG/CGF/Laboratory/LIMS/drop-box-prod/gwas_primaryqc/",
        descrition="For ``CGEMS/CCAD`` use only, if lims_upload is set to true, the LimsUpload file will be copied to LIMS production directory or designated directory.",
    )

    case_control_gwas: bool = Field(
        False,
        description="A plink logistic regression gwas will be performed with case_control phenotype.",
    )

    max_time_hr: Optional[int] = Field(
        None,
        description="The maximum amount of time that can be requested, in hours.",
    )

    max_mem_mb: Optional[int] = Field(
        None,
        description="The maximum amount of memory that can be requests, in megabytes.",
    )

    time_start: str = Field(
        timestr,
        description="Date and time at which the workflow starts. This creates a unique id for the run.",
    )
    convert_gtc2bcf: bool = Field(
        False,
        description="If input is GTC, this switches between gtc2vcf (True) and gtc2ped (False - default) for conversion to BED",
    )

    additional_params_for_gtc2bcf: str = Field(
        "--use-gtc-sample-names",
        description="Additional/optional parameters not hardcoded to be used or skipped in gtc2bcf for specific analysis",
    )

    @staticmethod
    def schema_rst():
        """Tweak schema for rendering in Sphinx."""
        import copy
        import json

        content = copy.deepcopy(WorkflowParams.schema())
        content["title"] = "Workflow Parameters"

        return json.dumps(content, indent=2)
