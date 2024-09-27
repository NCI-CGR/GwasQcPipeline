#!/usr/bin/env bash

ml  plink/1.9

plink \
--bfile subject_level/subjects \
--pheno subject_level/gwas.txt \
--ci 0.95 \
--assoc \
--out delivery/gwas_subjects \
--memory 20000 \
--threads 2 \
--allow-no-sex

ml -plink/1.9
