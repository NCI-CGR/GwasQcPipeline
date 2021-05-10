from cgr_gwas_qc import load_config

cfg = load_config()


rule extract_fingerprint_snps:
    """Extract fingerprint SNPs from a PLINK data and convert to a GRAF."""
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
    output:
        "{prefix}.fpg",
    log:
        "{prefix}.fpg.log",
    conda:
        cfg.conda("graf")
    shell:
        # GRAF returns an exit code of 1, this captures it so snakemake will actually run.
        "graf "
        "-exfp {wildcards.prefix} "
        "-out {output[0]} "
        "-type 4 "
        "> {log} 2>&1 "
        "|| exit_code=$?; if [ $exit_code -ne 1 ]; then exit $exit_code; fi"


rule relatedness:
    """Estimate relatedness using GRAF

    Outputs a table with pairwise samples and their genotypic relationship.
    GRAF uses a set of 10K pre-screened SNPs, so we can directly use Call
    Rate 2 filtered samples.
    """
    input:
        "{prefix}.fpg",
    output:
        "{prefix}_graf_relatedness.tsv",
    conda:
        cfg.conda("graf")
    log:
        "{prefix}_graf_relatedness.log",
    shell:
        # GRAF returns an exit code of 1, this captures it so snakemake will actually run.
        "graf "
        "-geno {input[0]} "
        "-type 4 "
        "-out {output[0]} "
        "> {log} 2>&1 "
        "|| exit_code=$?; if [ $exit_code -ne 1 ]; then exit $exit_code; fi"


rule relatedness_png:
    input:
        "sample_level/concordance/graf_relatedness.txt",
    output:
        "sample_level/plots/graf_relatedness.png",
    conda:
        cfg.conda("graf")
    shell:
        "PlotGraf.pl {input[0]} {output[0]} 3"


rule populations:
    """Estimate ancestry for each sample."""
    input:
        "{prefix}.fpg",
    output:
        "{prefix}_graf_populations.txt",
    log:
        "{prefix}_graf_populations.log",
    conda:
        cfg.conda("graf")
    shell:
        # GRAF returns an exit code of 1, this captures it so snakemake will actually run.
        "graf "
        "-geno {input[0]} "
        "-pop {output[0]} "
        "> {log} 2>&1 "
        "|| exit_code=$?; if [ $exit_code -ne 1 ]; then exit $exit_code; fi"


rule ancestry:
    """Create summary table with ancestry calls."""
    input:
        "{prefix}_graf_populations.txt",
    output:
        "{prefix}_graf_ancestry.txt",
    conda:
        cfg.conda("graf")
    shell:
        "PlotPopulations.pl {input[0]} {output[0]} "


rule populations_png:
    """Create summary table with ancestry calls."""
    input:
        "{prefix}_graf_populations.txt",
    output:
        "{prefix}_graf_populations.png",
    conda:
        cfg.conda("graf")
    shell:
        "PlotPopulations.pl {input[0]} {output[0]} "
