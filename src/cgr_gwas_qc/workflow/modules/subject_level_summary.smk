import pandas as pd


rule subject_qc_table:
    input:
        sample_qc_csv="sample_level/sample_qc.csv",
        sample_concordance_csv="sample_level/concordance/summary.csv",
    params:
        remove_sex_discordant=cfg.config.workflow_params.remove_sex_discordant,
        remove_unexpected_rep=cfg.config.workflow_params.remove_unexpected_rep,
    output:
        "subject_level/subject_qc.csv",
    script:
        "../scripts/subject_qc_table.py"


rule subject_representative:
    """Select samples to be subject representative.

    Pull out samples that have no analytic exclusions. Create mapping file
    for plink to convert Sample_IDs to Subject_IDs.
    """
    input:
        rules.subject_qc_table.output[0],
    output:
        sample2subject=temp("subject_level/sample_to_subject.txt"),
        selected=temp("subject_level/selected_subjects.txt"),
    run:
        from cgr_gwas_qc.workflow.scripts import subject_qc_table

        df = subject_qc_table.read(input[0]).query("not subject_analytic_exclusion")

        df.reindex(
            ["Sample_ID", "Sample_ID", "Group_By_Subject_ID", "Group_By_Subject_ID"], axis=1,
        ).to_csv(output.sample2subject, sep=" ", index=False, header=False)

        df.reindex(["Group_By_Subject_ID", "Group_By_Subject_ID"], axis=1).to_csv(
            output.selected, sep=" ", index=False, header=False
        )


rule convert_samples_to_subjects:
    input:
        bed="sample_level/call_rate_2/samples.bed",
        bim="sample_level/call_rate_2/samples.bim",
        fam="sample_level/call_rate_2/samples.fam",
        sample2subject=rules.subject_representative.output.sample2subject,
        selected=rules.subject_representative.output.selected,
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
        "--update-ids {input.sample2subject} "
        "--keep {input.selected} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


rule related_subjects:
    input:
        "{{prefix}}/subjects_maf{maf}_ld{ld}_pruned.concordance.csv".format(
            maf=cfg.config.software_params.maf_for_ibd, ld=cfg.config.software_params.ld_prune_r2,
        ),
    output:
        relatives="{prefix}/relatives.csv",
        to_remove="{prefix}/related_subjects_to_remove.txt",
    script:
        "../scripts/related_subjects.py"


rule remove_related_subjects:
    input:
        bed="{prefix}/subjects.bed",
        bim="{prefix}/subjects.bim",
        fam="{prefix}/subjects.fam",
        to_remove=rules.related_subjects.output.to_remove,
    params:
        out_prefix="{prefix}/subjects_unrelated",
    output:
        bed="{prefix}/subjects_unrelated.bed",
        bim="{prefix}/subjects_unrelated.bim",
        fam="{prefix}/subjects_unrelated.fam",
        nosex="{prefix}/subjects_unrelated.nosex",
    log:
        "{prefix}/subjects_unrelated.log",
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
