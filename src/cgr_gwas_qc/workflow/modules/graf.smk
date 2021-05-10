rule graf_extract_fingerprint_snps:
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
        "graf "
        "-exfp {wildcards.prefix} "
        "-out {output[0]} "
        "-type 4 "
        "> {log} 2>&1 "
        "|| exit_code=$?; if [ $exit_code -ne 1 ]; then exit $exit_code; fi" # GRAF returns an exit code of 1, this captures it so snakemake will actually run.


rule graf_relatedness_png:
    input:
        "sample_level/concordance/graf_relatedness.txt",
    output:
        "sample_level/plots/graf_relatedness.png",
    conda:
        cfg.conda("graf")
    shell:
        "PlotGraf.pl {input[0]} {output[0]} 3"


rule graf_ancestry_png:
    """Create summary table with ancestry calls."""
    input:
        "sample_level/ancestry/graf_populations.txt",
    output:
        "sample_level/plots/graf_ancestry.png",
    conda:
        cfg.conda("graf")
    shell:
        "PlotPopulations.pl {input[0]} {output[0]} "
