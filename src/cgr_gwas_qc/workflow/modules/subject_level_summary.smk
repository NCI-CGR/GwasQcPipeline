import pandas as pd


rule samples_to_keep:
    input:
        "sample_level/samples_used_for_subjects.csv",
    output:
        temp("subject_level/samples_to_keep.txt"),
    run:
        (
            pd.read_csv(input[0])
            .dropna()
            .assign(Sample_ID2=lambda x: x.Sample_ID)
            .reindex(["Sample_ID", "Sample_ID2"], axis=1)
            .to_csv(output[0], sep=" ", index=False, header=False)
        )


rule kept_samples:
    input:
        bed="sample_level/call_rate_2/samples.bed",
        bim="sample_level/call_rate_2/samples.bim",
        fam="sample_level/call_rate_2/samples.fam",
        to_keep=rules.samples_to_keep.output[0],
    params:
        out_prefix="subject_level/samples",
    output:
        bed=temp("subject_level/samples.bed"),
        bim=temp("subject_level/samples.bim"),
        fam=temp("subject_level/samples.fam"),
    log:
        "subject_level/samples.log",
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


rule sample_to_subject_map:
    input:
        "sample_level/samples_used_for_subjects.csv",
    output:
        temp("subject_level/samples_to_subjects.txt"),
    run:
        (
            pd.read_csv(input[0])
            .dropna()
            .assign(Sample_ID2=lambda x: x.Sample_ID)
            .assign(Subject_ID2=lambda x: x.Subject_ID)
            .reindex(["Sample_ID", "Sample_ID2", "Subject_ID", "Subject_ID2"], axis=1)
            .to_csv(output[0], sep=" ", index=False, header=False)
        )


rule renamed_subjects:
    input:
        bed=rules.kept_samples.output.bed,
        bim=rules.kept_samples.output.bim,
        fam=rules.kept_samples.output.fam,
        to_rename=rules.sample_to_subject_map.output[0],
    params:
        out_prefix="subject_level/subjects",
    output:
        bed="subject_level/subjects.bed",
        bim="subject_level/subjects.bim",
        fam="subject_level/subjects.fam",
    log:
        "subject_level/subjects.log",
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
        "--update-ids {input.to_rename} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem} "
        "--out {params.out_prefix}"
