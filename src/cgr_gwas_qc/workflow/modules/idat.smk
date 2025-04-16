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
            red = cfg.expand(
                cfg.config.user_files.idat_pattern.red, query="is_missing_idats==False"
            )
            green = cfg.expand(
                cfg.config.user_files.idat_pattern.green, query="is_missing_idats==False"
            )
        else:
            red = cfg.expand(
                cfg.config.user_files.idat_pattern.red,
                query=f'cluster_group=="{wildcards.grp}"&is_missing_idats==False',
            )
            green = cfg.expand(
                cfg.config.user_files.idat_pattern.green,
                query=f'cluster_group=="{wildcards.grp}"&is_missing_idats==False',
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
    threads: workflow.cores
    envmodules:
        "dragena/1.0.0",
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * 100 * attempt,
        time_hr=lambda wildcards, attempt: 5 * attempt,
    shell:
        """
        if [ "{params.dragena_location}" != "None" ];then dragena='{params.dragena_location}';else  dragena='dragena';fi
        if ! command -v $dragena &> /dev/null; then
            echo "Error: Dragena executable not found. Ensure that 'dragena/1.0.0' module is loaded or provide a valid 'dragena_location' in the config." >&2
            exit 1
        else
            $dragena genotype call --bpm-manifest {params.bpm} --cluster-file {params.egt} --idat-sample-sheet {input.idat_ss} --num-threads {threads} --output-folder {output.output_folder}
        fi
        """


rule check_gtc_creation:
    """Checks if GTC files are created. If not, it will mark the sample as missing GTC.
    """
    input:
        "sample_level/gtcs",
    output:
        "sample_level/gtcs_check.done",
    run:
        from cgr_gwas_qc import load_config
        import pandas as pd
        from pathlib import Path
        from itertools import chain

        cfg = load_config()


        def get_gtcs(input):
            gtcList = [str(p.stem) for p in Path(input).rglob("*.gtc")]
            return gtcList


        if isinstance(input, list):
            gtcList = list(chain(*(map(get_gtcs, input))))
        else:
            gtcList = get_gtcs(input[0])
        cfg.ss.is_missing_gtc = pd.Series(
            barcode not in gtcList
            for barcode in cfg.expand("{SentrixBarcode_A}_{SentrixPosition_A}")
        )
        cfg.ss.to_csv("cgr_sample_sheet.csv", index=False)
        Path("sample_level/gtcs_check.done").touch()
