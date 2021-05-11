from cgr_gwas_qc import load_config

cfg = load_config()


rule convert:
    input:
        ped="{prefix}.ped",
        map_="{prefix}.map",
    output:
        par=temp("{prefix}.convert.par"),
        gen=temp("{prefix}.gen"),
        snp=temp("{prefix}.snp"),
        ind=temp("{prefix}.ind"),
    conda:
        cfg.conda("eigensoft")
    shell:
        """
        echo 'genotypename: {input.ped}' > {output.par} \
        && echo 'snpname: {input.map_}' >> {output.par} \
        && echo 'indivname: {input.ped}' >> {output.par} \
        && echo 'outputformat: EIGENSTRAT' >> {output.par} \
        && echo 'genooutfilename: {output.gen}' >> {output.par} \
        && echo 'snpoutfilename: {output.snp}' >> {output.par} \
        && echo 'indoutfilename: {output.ind}' >> {output.par} \
        && echo 'familynames: NO' >> {output.par} \
        && convertf -p {output.par}
        """


rule smartpca:
    input:
        gen="{prefix}.gen",
        snp="{prefix}.snp",
        ind="{prefix}.ind",
    output:
        par=temp("{prefix}.pca.par"),
        eigenvec="{prefix}.eigenvec",
    conda:
        cfg.conda("eigensoft")
    shell:
        """
        echo 'genotypename: {input.gen}' > {output.par} \
        && echo 'snpname: {input.snp}' >> {output.par} \
        && echo 'indivname: {input.ind}' >> {output.par} \
        && echo 'evecoutname: {output.eigenvec}' >> {output.par} \
        && echo 'fastmode: YES >> {output.par}' \
        && smartpca -p {output.par}
        """
