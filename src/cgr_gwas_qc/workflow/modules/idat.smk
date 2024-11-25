from cgr_gwas_qc import load_config

cfg = load_config()


rule write_idat2gtc_ss:
    """Writes a samplesheet for Dragen array idat to gtc.
    """
    params:
        grp="",
    output:
        "sample_level/idat.csv",
    run:
        import pandas as pd

        if params.grp == "":
            red = cfg.expand(cfg.config.user_files.idat_pattern.red)
            green = cfg.expand(cfg.config.user_files.idat_pattern.green)
        else:
            red = cfg.expand(
                cfg.config.user_files.idat_pattern.red, query=f'cluster_group=="{wildcards.grp}"'
            )
            green = cfg.expand(
                cfg.config.user_files.idat_pattern.green, query=f'cluster_group=="{wildcards.grp}"'
            )
        pd.DataFrame({"Green IDAT Path": green, "Red IDAT Path": red}).to_csv(
            output[0], index=False
        )


rule idat2gtc:
    input:
        idat_ss=rules.write_idat2gtc_ss.output[0],
    params:
        bpm=cfg.config.reference_files.illumina_manifest_file,
        egt=cfg.config.reference_files.illumina_cluster_file,
        dragena_location=cfg.config.workflow_params.dragena_location,
    output:
        output_folder=directory("sample_level/gtcs"),
    threads: 44
    envmodules:
        "dragena/1.0.0",
    shell:
        """
        if [ "{params.dragena_location}" != "None" ];then dragena='{params.dragena_location}';else  dragena='dragena';fi
        $dragena genotype call --bpm-manifest {params.bpm} --cluster-file {params.egt} --idat-sample-sheet {input.idat_ss} --num-threads {threads} --output-folder {output.output_folder}
        """
