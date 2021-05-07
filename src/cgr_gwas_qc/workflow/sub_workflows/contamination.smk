"""Sample Contamination Workflow

Inputs:
    - IDAT Files (red and green)
    - GTC Files
    - Illumina BPM
    - 1000 Genomes VCF

Outputs:

    - sample_level/contamination.csv
"""
from cgr_gwas_qc import load_config

cfg = load_config()


include: cfg.modules("plink_filters")
include: cfg.modules("plink_stats")


def _contamination_outputs(wildcards):
    if (
        cfg.config.user_files.idat_pattern
        and cfg.config.user_files.gtc_pattern
        and cfg.config.workflow_params.remove_contam
    ):
        return [
            "sample_level/contamination.csv",
            "sample_level/median_idat_intensity.csv",
        ]

    return []


rule all_contamination:
    input:
        _contamination_outputs,


rule per_sample_median_idat_intensity:
    """Calculate median intensity of Red and Green channels."""
    input:
        red=lambda wc: cfg.expand(
            cfg.config.user_files.idat_pattern.red, query=f"Sample_ID == '{wc.Sample_ID}'",
        ),
        green=lambda wc: cfg.expand(
            cfg.config.user_files.idat_pattern.green, query=f"Sample_ID == '{wc.Sample_ID}'",
        ),
    output:
        temp(
            "sample_level/per_sample_median_idat_intensity/{Sample_ID}.{SentrixBarcode_A}.{SentrixPosition_A}.csv"
        ),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
    group:
        "per_sample_median_idat_intensity"
    conda:
        cfg.conda("illuminaio")
    script:
        "../scripts/median_idat_intensity.R"


rule agg_median_idat_intensity:
    input:
        cfg.expand(rules.per_sample_median_idat_intensity.output[0]),
    output:
        "sample_level/median_idat_intensity.csv",
    run:
        pd.concat([pd.read_csv(file_name) for file_name in input]).to_csv(output[0], index=False)


rule per_sample_gtc_to_adpc:
    """Converts a sample's GTC/BPM to an Illumina ADPC.BIN.

    This is the format required by ``verifyIDintensity``. The script also
    runs some sanity checks (intensities and normalized intensities > 0;
    genotypes are one of {0, 1, 2, 3}) while processing each file.

    .. warning::
        This is a submission hot spot creating 1 job per sample.
    """
    input:
        gtc_file=lambda wc: cfg.expand(
            cfg.config.user_files.gtc_pattern, query=f"Sample_ID == '{wc.Sample_ID}'",
        )[0],
        bpm_file=cfg.config.reference_files.illumina_manifest_file,
    output:
        temp("sample_level/per_sample_adpc/{Sample_ID}.adpc.bin"),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
    group:
        "per_sample_gtc_to_adpc"
    script:
        "../scripts/gtc2adpc.py"


rule pull_1KG_allele_b_freq:
    """Pulls the population level allele frequencies from the 1KG project.

    ``verifyIDintensity`` requires population level allele frequencies
    for its model. Here we use a custom script to pull out the allele B
    frequencies (ABF) from the 1000 genomes project (1KG). To do this we
    take each marker from the manifest file (BPM) and pull out ABF in the
    1KG ``.vcf`` from the ``INFO`` column. The script allows pulling out
    allele frequencies for different super populations but defaults to
    ``AF`` which ignores super population.
    """
    input:
        bpm_file=cfg.config.reference_files.illumina_manifest_file,
        vcf_file=cfg.config.reference_files.thousand_genome_vcf,
    params:
        population=cfg.config.software_params.contam_population,
    output:
        abf_file="sample_level/{}.{}.abf.txt".format(
            cfg.config.reference_files.illumina_manifest_file.stem,
            cfg.config.software_params.contam_population,
        ),
    script:
        "../scripts/bpm2abf.py"


rule per_sample_verifyIDintensity_contamination:
    """Find contaminated samples using allele intensities.

    Uses ``verifyIDintensity`` to find samples with allele intensities that deviate from the
    population.

    .. warning::
        This is a submission hot spot creating 1 job per sample.

    .. note::
        Here we are running ``verifyIDintensity`` in single sample mode. This software also has a
        multi-sample mode which may be faster and give better estimates. The problem with
        multi-sample mode is that it only works when you have a "large" number of samples.
    """
    input:
        adpc=rules.per_sample_gtc_to_adpc.output[0],
        abf=rules.pull_1KG_allele_b_freq.output.abf_file,
    params:
        snps=cfg.config.num_snps,
    output:
        temp("sample_level/per_sample_contamination_test/{Sample_ID}.contam.out"),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
    group:
        "per_sample_verifyIDintensity_contamination"
    conda:
        cfg.conda("verifyidintensity")
    shell:
        "verifyIDintensity -m {params.snps} -n 1 -b {input.abf} -v -p -i {input.adpc} > {output}"


rule agg_verifyIDintensity_contamination:
    """Aggregate sample contamination scores.

    Aggregates sample level contamination scores into a single file (each
    row is a sample). The script sets ``%Mix`` to ``NA`` if the intensity
    is below the threshold and the file is not in the ``imiss`` file.
    """
    input:
        contamination_files=cfg.expand(rules.per_sample_verifyIDintensity_contamination.output),
        median_intensity_file=rules.agg_median_idat_intensity.output[0],
        imiss_file="sample_level/call_rate_2/samples.imiss",
    params:
        intensity_threshold=cfg.config.software_params.intensity_threshold,
        contam_threshold=cfg.config.software_params.contam_threshold,
    output:
        "sample_level/contamination.csv",
    script:
        "../scripts/agg_contamination.py"
