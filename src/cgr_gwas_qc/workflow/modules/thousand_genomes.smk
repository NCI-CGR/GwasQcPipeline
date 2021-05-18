from cgr_gwas_qc import load_config

cfg = load_config()


rule pull_b_allele_freq_from_1kg:
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


rule update_snps_to_1kg_rsID:
    """Update SNP IDs to rsID from the 1KG project.

    Update study marker IDs to correspond with the 1KG project rsIDs.
    """
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
        vcf=cfg.config.reference_files.thousand_genome_vcf,
    output:
        bed="{prefix}_1kg_rsID.bed",
        bim="{prefix}_1kg_rsID.bim",
        fam="{prefix}_1kg_rsID.fam",
        id_map="{prefix}_1kg_rsID.csv",
    script:
        "../scripts/update_snps_to_1kg_rsID.py"
