from pathlib import Path
from textwrap import dedent

import pandas as pd


wildcard_constraints:
    name="samples|subjects|controls|snps",
    filters=".*",
    cr="1|2",
    maf="0\.\d+",
    ld="0\.\d+",
    deliver_prefix=".*",
    deliver_suffix=".*",


def rule_group(wildcards):
    prefix = wildcards.get("prefix", "").replace("/", "_")
    name = wildcards.get("name", "")
    return f"{prefix}_{name}"


rule concordance_table:
    """Parse IBD file and calculate sample concordance.

    Returns:
        DataFrame with ["ID1", "ID2", "PI_HAT", "concordance"].
          Sample concordance is `IBS2 / (IBS0 + IBS1 + IBS2)`.
    """
    input:
        "{prefix}.genome",
    params:
        concordance_threshold=cfg.config.software_params.dup_concordance_cutoff,
        pi_hat_threshold=cfg.config.software_params.pi_hat_threshold,
    output:
        "{prefix}.concordance.csv",
    script:
        "../scripts/concordance_table.py"
