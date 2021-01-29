
rule graf_relatedness_png:
    input:
        "sample_level/concordance/graf_relatedness.txt",
    output:
        "sample_level/plots/graf_relatedness.png",
    conda:
        cfg.conda("graf.yml")
    shell:
        "PlotGraf.pl {input[0]} {output[0]} 3"


rule graf_ancestry_png:
    """Create summary table with ancestry calls."""
    input:
        "sample_level/ancestry/graf_populations.txt",
    output:
        "sample_level/plots/graf_ancestry.png",
    conda:
        cfg.conda("graf.yml")
    shell:
        "PlotPopulations.pl {input[0]} {output[0]} "
