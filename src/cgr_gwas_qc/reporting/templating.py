from jinja2 import Environment, PackageLoader

env = Environment(
    loader=PackageLoader("cgr_gwas_qc", "reporting/templates"),
    trim_blocks=True,
    lstrip_blocks=True,
    keep_trailing_newline=True,
)


def number_formater(value, decimals: int = 2) -> str:
    if isinstance(value, float):
        return f"{value:,.{decimals}f}"

    if isinstance(value, int):
        return f"{value:,}"

    return str(value)


def to_pct(value) -> float:
    return value * 100


env.filters["numFormat"] = number_formater
env.filters["toPct"] = to_pct
