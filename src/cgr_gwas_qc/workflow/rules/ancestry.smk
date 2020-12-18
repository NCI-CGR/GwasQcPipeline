from cgr_gwas_qc.testing.data import Graf

graf = Graf()


rule graf_fingerprint_list:
    input:
        graf / "data/G1000FpGeno.bim",
    output:
        "ancestry/graf_fingerprints.txt",
    shell:
        "awk '{{print $2}}' {input[0]} > {output[0]}"


rule rename_to_thousG:
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


rule convert_to_fpg:
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
    input:
        rules.convert_to_fpg.output[0],
    output:
        "ancestry/graf_relatedness.txt",
    shell:
        f"export PATH={str(graf)}:$PATH; "
        "graf "
        "-geno {input[0]} "
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
    input:
        rules.convert_to_fpg.output[0],
    output:
        "ancestry/graf_pop.txt",
    shell:
        f"export PATH={str(graf)}:$PATH; "
        "graf "
        "-geno {input[0]} "
        "-pop {output[0]} "

        "|| exit_code=$?; if [ $exit_code -ne 1 ]; then exit $exit_code; fi" # GRAF returns an exit code of 1, this captures it so snakemake will actually run.


rule graf_make_ancestry_call_table:
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
