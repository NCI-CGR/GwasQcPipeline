import pandas as pd


rule subject_representative:
    input:
        "sample_level/sample_qc.csv",
    output:
        temp("subject_level/subject_representative.txt"),
    run:
        (
            pd.read_csv(input[0])
            .query("is_subject_representative")
            .assign(Sample_ID2=lambda x: x.Sample_ID)
            .reindex(["Sample_ID", "Sample_ID2"], axis=1)
            .to_csv(output[0], sep=" ", index=False, header=False)
        )


rule kept_samples:
    input:
        bed="sample_level/call_rate_2/samples.bed",
        bim="sample_level/call_rate_2/samples.bim",
        fam="sample_level/call_rate_2/samples.fam",
        to_keep=rules.subject_representative.output[0],
    params:
        out_prefix="subject_level/samples",
    output:
        bed=temp("subject_level/samples.bed"),
        bim=temp("subject_level/samples.bim"),
        fam=temp("subject_level/samples.fam"),
        nosex=temp("subject_level/samples.nosex"),
    log:
        "subject_level/samples.log",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--keep {input.to_keep} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


rule sample_to_subject_map:
    input:
        "sample_level/sample_qc.csv",
    output:
        temp("subject_level/samples_to_subjects.txt"),
    run:
        (
            pd.read_csv(input[0])
            .query("is_subject_representative")
            .assign(Sample_ID2=lambda x: x.Sample_ID)
            .assign(Subject_ID2=lambda x: x.Group_By_Subject_ID)
            .reindex(["Sample_ID", "Sample_ID2", "Group_By_Subject_ID", "Subject_ID2"], axis=1)
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
        nosex="subject_level/subjects.nosex",
    log:
        "subject_level/subjects.log",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--update-ids {input.to_rename} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


rule related_subjects:
    input:
        ibs="{{prefix}}/subjects_maf{maf}_ld{ld}_pruned.genome".format(
            maf=cfg.config.software_params.maf_for_ibd, ld=cfg.config.software_params.ld_prune_r2,
        ),
    params:
        pi_hat_threshold=lambda wc: float(wc.pi),
    output:
        "{prefix}/subjects_to_remove_pi_hat_gt{pi}.txt",
    wildcard_constraints:
        pi="[01].\d+",
    script:
        "../scripts/related_subjects.py"


rule remove_related_subjects:
    input:
        bed="{prefix}/subjects.bed",
        bim="{prefix}/subjects.bim",
        fam="{prefix}/subjects.fam",
        to_remove=rules.related_subjects.output[0],
    params:
        out_prefix="{prefix}/subjects_unrelated{pi}",
    output:
        bed="{prefix}/subjects_unrelated{pi}.bed",
        bim="{prefix}/subjects_unrelated{pi}.bim",
        fam="{prefix}/subjects_unrelated{pi}.fam",
        nosex="{prefix}/subjects_unrelated{pi}.nosex",
    log:
        "{prefix}/subjects_unrelated{pi}.log",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--remove {input.to_remove} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"
