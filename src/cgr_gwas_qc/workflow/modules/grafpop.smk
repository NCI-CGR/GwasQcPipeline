from cgr_gwas_qc import load_config

cfg = load_config()

rule extract_fingerprint_snps:
    """Extract fingerprint SNPs from a PLINK data and convert to a GRAF."""
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
    params:
        out_prefix="{prefix}",
    output:
        "{prefix}.fpg",
    log:
        "{prefix}.fpg.log",
    conda:
        cfg.conda("graf")
    shell:
        # GRAF returns an exit code of 1, this captures it so snakemake will actually run.
        "graf "
        "-exfp {params.out_prefix} "
        "-out {output[0]} "
        "-type 4 "
        "> {log} 2>&1 "
        "|| exit_code=$?; if [ $exit_code -ne 1 ]; then exit $exit_code; fi"


rule grafpop_populations:
    """Estimate ancestry for each sample using grafpop1.0."""
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
    output:
        "{prefix}_graf_populations.txt",
    log:
        "{prefix}_graf_populations.log",
    conda:
        cfg.conda("grafpop")
    shell:
        """
        # > in grafpop needed b/c graf has exit status 1 which causes snakemake to error out
        export GRAFPATH=$CONDA_PREFIX/share/grafpop
        grafpop {input.bim} {output[0]} > {log} 2>&1 || exit_code=$?; if [ $exit_code -ne 1 ]; then exit $exit_code; fi
        """

