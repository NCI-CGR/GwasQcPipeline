def schema_to_dict(schema: dict) -> dict:
    """Parse a pydantic schema into an dict with all default values."""

    def parse_ref(reference: dict) -> dict:
        ref = reference["$ref"].split("/")[-1]
        return schema["definitions"][ref]["properties"]

    def parse_properties(properties: dict, comments: bool = False) -> dict:
        res = {}
        for k, v in properties.items():
            if "$ref" in v.keys():
                v = parse_ref(v)
                res[k] = parse_properties(v)
            else:
                res[k] = v.get("description", "") if comments else v.get("default", "")
        return res

    return parse_properties(schema["properties"])
