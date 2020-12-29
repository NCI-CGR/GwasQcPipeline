from cgr_gwas_qc.testing.data import Graf

graf = Graf()


rule graf_fingerprint_list:
    """Create a list of GRAF's fingerprint markers."""
    input:
        graf / "data/G1000FpGeno.bim",
    output:
        "ancestry/graf_fingerprints.txt",
    shell:
        "awk '{{print $2}}' {input[0]} > {output[0]}"


rule rename_to_thousG:
    """Update SNP IDs to rsID from the 1KG project.

    Update study marker IDs to correspond with GRAF's fingerprints which are
    based on the 1KG project rsIDs.
    """
    input:
        bim="plink_filter_call_rate_2/samples.bim",
        vcf=cfg.config.reference_files.thousand_genome_vcf,
    output:
        snps_to_remove="ancestry/snps_not_in_thousG.txt",
        bim="ancestry/samples.renamed_to_thousG.bim",
    log:
        "ancestry/rename_to_thousG.log",
    script:
        "../scripts/bim_filter_vcf.py"


rule extract_graf_fingerprint_markers:
    """Use PLINK to pull out GRAF's fingerprints rsIDs from study data."""
    input:
        bed="plink_filter_call_rate_2/samples.bed",
        bim=rules.rename_to_thousG.output.bim,
        fam="plink_filter_call_rate_2/samples.fam",
        graf_fp=rules.graf_fingerprint_list.output[0],
    params:
        in_prefix="plink_filter_call_rate_2/samples",
        out_prefix="ancestry/samples",
    output:
        bed="ancestry/samples.bed",
        bim="ancestry/samples.bim",
        fam="ancestry/samples.fam",
    log:
        "ancestry/samples.log",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    threads: 2
    resources:
        mem=10000,
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--extract {input.graf_fp} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem} "
        "--out {params.out_prefix}"


rule create_ssm:
    """Create mapping of Subject_ID to Sample_ID.

    GRAF has an option to pass this file, but I don't see where it is making
    any changes to the outputs. I am keeping it here because I do use this
    file during testing.
    """
    output:
        "ancestry/ssm.txt",
    run:
        (
            cfg.ss[["Group_By_Subject_ID", "Sample_ID"]]
            .rename({"Group_By_Subject_ID": "Subject_ID"}, axis=1)
            .dropna()
            .to_csv(output[0], sep="\t", index=False)
        )


rule convert_to_fpg:
    """Extract PLINK data set to GRAF data set."""
    input:
        bed=rules.extract_graf_fingerprint_markers.output.bed,
        bim=rules.extract_graf_fingerprint_markers.output.bim,
        fam=rules.extract_graf_fingerprint_markers.output.fam,
    params:
        in_prefix="ancestry/samples",
    output:
        "ancestry/samples.fpg",
    shell:
        f"export PATH={str(graf)}:$PATH; "
        "graf "
        "-exfp {params.in_prefix} "
        "-out {output[0]} "
        "-type 4 "
        "|| exit_code=$?; if [ $exit_code -ne 1 ]; then exit $exit_code; fi" # GRAF returns an exit code of 1, this captures it so snakemake will actually run.


rule graf_relatedness:
    """Estimate relatedness among samples.

    Outputs a table with pairwise samples and their genotypic relationship.
    """
    input:
        fpg=rules.convert_to_fpg.output[0],
        ssm=rules.create_ssm.output[0],
    output:
        "ancestry/graf_relatedness.txt",
    shell:
        f"export PATH={str(graf)}:$PATH; "
        "graf "
        "-geno {input.fpg} "
        "-ssm {input.ssm} "
        "-out {output[0]} "
        "|| exit_code=$?; if [ $exit_code -ne 1 ]; then exit $exit_code; fi" # GRAF returns an exit code of 1, this captures it so snakemake will actually run.


rule graf_relatedness_png:
    input:
        rules.graf_relatedness.output[0],
    params:
        script=(graf / "PlotGraf.pl").as_posix(),
    output:
        "ancestry/graf_relatedness.png",
    conda:
        cfg.conda("graf_perl.yml")
    shell:
        "perl {params.script} {input[0]} {output[0]} 3"


rule graf_ancestry:
    """Estimate ancestry for each sample."""
    input:
        fpg=rules.convert_to_fpg.output[0],
        ssm=rules.create_ssm.output[0],
    output:
        "ancestry/graf_pop.txt",
    shell:
        f"export PATH={str(graf)}:$PATH; "
        "graf "
        "-geno {input.fpg} "
        "-ssm {input.ssm} "
        "-pop {output[0]} "
        "|| exit_code=$?; if [ $exit_code -ne 1 ]; then exit $exit_code; fi" # GRAF returns an exit code of 1, this captures it so snakemake will actually run.


rule graf_make_ancestry_call_table:
    """Create summary table with ancestry calls."""
    input:
        rules.graf_ancestry.output[0],
    params:
        script=(graf / "PlotPopulations.pl").as_posix(),
    output:
        calls="ancestry/graf_ancestry_calls.txt",
        png="ancestry/graf_ancestry_calls.png",
    conda:
        cfg.conda("graf_perl.yml")
    shell:
        "perl {params.script} {input[0]} {output.calls} "
        "&& perl {params.script} {input[0]} {output.png} "
