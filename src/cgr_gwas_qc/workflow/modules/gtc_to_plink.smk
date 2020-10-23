"""Convert genotype calls to a ``plink`` compatible format.

Converts GTC + BPM for all samples and then merges them into the binary
``plink`` format (BED + BIM + FAM).
"""


rule gtc_to_ped:
    """Converts a sample's GTC + BPM to PED + MAP.

    .. warning::
        This is a submission hot spot creating 1 job per sample.
    """
    input:
        gtc=lambda wc: cfg.expand(
            cfg.config.file_patterns.gtc, query=f"Sample_ID == '{wc.Sample_ID}'"
        )[0],
        bpm=cfg.config.reference_files.illumina_manifest_file,
    params:
        strand=cfg.config.software_params.strand,
    output:
        ped=temp("per_sample_plink_files/{Sample_ID}.ped"),
        map_=temp("per_sample_plink_files/{Sample_ID}.map"),
    script:
        "../scripts/gtc2plink.py"


rule sample_ped_merge_list:
    input:
        ped=cfg.expand("per_sample_plink_files/{Sample_ID}.ped"),
        map_=cfg.expand("per_sample_plink_files/{Sample_ID}.map"),
    output:
        "plink_start/mergeList.txt",
    group:
        "merge_sample_peds"
    run:
        with open(output[0], "w") as fh:
            for ped, map_ in zip(input.ped, input.map_):
                fh.write(f"{ped} {map_}\n")


rule merge_sample_peds:
    """Merges multiple ``plink`` samples.

    Merge and convert sample(s) into the ``plink`` binary format.
    """
    input:
        ped=cfg.expand("per_sample_plink_files/{Sample_ID}.ped"),
        map_=cfg.expand("per_sample_plink_files/{Sample_ID}.map"),
        merge_list=rules.sample_ped_merge_list.output[0],
    params:
        prefix="plink_start/samples",
    output:
        bed="plink_start/samples.bed",
        bim="plink_start/samples.bim",
        fam="plink_start/samples.fam",
    log: "plink_start/samples.log",
    group:
        "merge_sample_peds"
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    threads: 20
    resources:
        mem=24000,
    shell:
        "plink --merge-list {input.merge_list} --make-bed --out {params.prefix} --threads {threads} --memory {resources.mem}"
