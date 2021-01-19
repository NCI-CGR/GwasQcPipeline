import pandas as pd


rule subjects_per_population:
    input:
        "sample_level/qc_summary.csv",
    output:
        "population_level/subject_list_{population}.txt",
    run:
        (
            pd.read_csv(input[0])
            .query("Subject_Representative & Ancestry == @wildcards.population")
            .assign(Subject_ID2=lambda x: x.Group_By_Subject_ID)
            .reindex(["Group_By_Subject_ID", "Subject_ID2"], axis=1)
            .to_csv(output[0], sep=" ", index=False)
        )


rule plink_split_population:
    input:
        bed="subject_level/subjects.bed",
        bim="subject_level/subjects.bim",
        fam="subject_level/subjects.fam",
        to_keep=rules.subjects_per_population.output[0],
    params:
        out_prefix="population_level/{population}",
    output:
        bed="population_level/{population}.bed",
        bim="population_level/{population}.bim",
        fam="population_level/{population}.fam",
        nosex="population_level/{population}.nosex",
    log:
        "population_level/{population}.log",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    threads: 20
    resources:
        mem=10000,
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--keep {input.to_keep} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem} "
        "--out {params.out_prefix}"
