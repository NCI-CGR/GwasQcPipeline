from cgr_gwas_qc import load_config

cfg = load_config()


def eigensoft_config_inputs(wildcards):
    tool = wildcards.tool
    prefix = wildcards.prefix
    if tool == "pca":
        return {
            "gen": f"{prefix}.gen",
            "snp": f"{prefix}.snp",
            "ind": f"{prefix}.ind",
        }

    # Default is to convert PED/MAP to EIGENSTRAT
    return {
        "gen": f"{prefix}.ped",
        "snp": f"{prefix}.map",
        "ind": f"{prefix}.ped",
    }


def eigensoft_config_params(wildcards):
    tool = wildcards.tool
    prefix = wildcards.prefix

    if tool == "pca":
        return f"""\
        genotypename: {prefix}.gen
        snpname: {prefix}.snp
        indivname: {prefix}.ind
        evecoutname: {prefix}.eigenvec
        fastmode: YES
        """

    # Default is to convert PED/MAP to EIGENSTRAT
    return f"""\
    genotypename: {prefix}.ped
    snpname: {prefix}.map
    indivname: {prefix}.ped
    outputformat: EIGENSTRAT
    genooutfilename: {prefix}.gen
    snpoutfilename: {prefix}.snp
    indoutfilename: {prefix}.ind
    familynames: NO
    """


rule eigensoft_create_par:
    input:
        unpack(eigensoft_config_inputs),
    params:
        eigensoft_config_params,
    output:
        temp("{prefix}.{tool}.par"),
    wildcard_constraints:
        tool="convert|pca",
    run:
        from textwrap import dedent

        Path(output[0]).write_text(dedent(params[0]))


rule eigensoft_convert:
    input:
        ped="{prefix}.ped",
        map_="{prefix}.map",
        par="{prefix}.convert.par",
    output:
        gen=temp("{prefix}.gen"),
        snp=temp("{prefix}.snp"),
        ind=temp("{prefix}.ind"),
    conda:
        cfg.conda("eigensoft")
    shell:
        "convertf -p {input.par}"


rule eigensoft_smartpca:
    input:
        gen="{prefix}.gen",
        snp="{prefix}.snp",
        ind="{prefix}.ind",
        par="{prefix}.pca.par",
    output:
        "{prefix}.eigenvec",
    conda:
        cfg.conda("eigensoft")
    shell:
        "smartpca -p {input.par}"
