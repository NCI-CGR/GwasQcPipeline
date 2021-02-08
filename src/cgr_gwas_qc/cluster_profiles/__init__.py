from jinja2 import Environment, PackageLoader

env = Environment(
    loader=PackageLoader("cgr_gwas_qc", "cluster_profiles"),
    trim_blocks=True,
    lstrip_blocks=True,
    keep_trailing_newline=True,
)
