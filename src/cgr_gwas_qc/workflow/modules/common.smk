from pathlib import Path
from textwrap import dedent

from cgr_gwas_qc.parsers.illumina import BeadPoolManifest


numSNPs = BeadPoolManifest(cfg.config.reference_files.illumina_manifest_file).num_loci


rule plink_bed_to_ped:
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
    output:
        ped="{prefix}.ped",
        map_="{prefix}.map",
    log:
        "{prefix}.log",
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
        "--recode "
        "--keep-allele-order "
        "--threads {threads} "
        "--memory {resources.mem} "
        "--out {wildcards.prefix}"


rule eigenstrat_config:
    input:
        ped="{prefix}.ped",
        map_="{prefix}.map",
    params:
        gen="{prefix}.eigenstratgeno",
        snp="{prefix}.snp",
        ind="{prefix}.ind",
    output:
        par="{prefix}.convertEigen.par",
    group:
        "eigenstrat"
    run:
        Path(output.par).write_text(
            dedent(
                f"""\
        genotypename: {input.ped}
        snpname: {input.map_}
        indivname: {input.ped}
        outputformat: EIGENSTRAT
        genooutfilename: {params.gen}
        snpoutfilename: {params.snp}
        indoutfilename: {params.ind}
        familynames: NO
        """
            )
        )


rule eigenstrat_convert:
    input:
        ped="{prefix}.ped",
        map_="{prefix}.map",
        par="{prefix}.convertEigen.par",
    output:
        gen="{prefix}.eigenstratgeno",
        snp="{prefix}.snp",
        ind="{prefix}.ind",
    group:
        "eigenstrat"
    envmodules:
        cfg.envmodules("eigensoft"),
    conda:
        cfg.conda("eigensoft.yml")
    shell:
        "convertf -p {input.par}"
