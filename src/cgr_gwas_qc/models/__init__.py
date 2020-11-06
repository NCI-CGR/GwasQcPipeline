def schema_to_dict(schema: dict) -> dict:
    """Parse a pydantic schema into an dict with all default values.

    Scans through a pydantic schema and creates a simplified dictionary where
    all values are replaced with their defaults or an empty string. This is
    useful for creating a ``config.yml``.

    Example:
        >>> from cgr_gwas_qc.models.config import Config
        >>> res = schema_to_dict(Config.schema())
        >>> res.reference_files.thousand_genome_vcf
        "/DCEG/CGF/Bioinformatics/Production/data/thousG/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz"
    """

    def parse_ref(reference: dict) -> dict:
        ref = reference["$ref"].split("/")[-1]
        return schema["definitions"][ref]["properties"]

    def parse_properties(properties: dict, comments: bool = False) -> dict:
        res = {}
        for k, v in properties.items():
            if "$ref" in v.keys():
                v = parse_ref(v)
                res[k] = parse_properties(v, comments)
            elif v.get("default", "") == "Deprecated":
                continue
            else:
                res[k] = v.get("description", "") if comments else v.get("default", "")
        return res

    return parse_properties(schema["properties"])
