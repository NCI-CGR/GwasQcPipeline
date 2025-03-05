# This script reads in eigensoft formatted snps and related samples
# smartsnp is used to calculate pca on unrelated subjects and project related subjects

require(smartsnp)
require(genio)
require(dplyr)

# read in data
geno_loc <- snakemake@input[['gen']]
related_samp <- read.table(snakemake@input[['related_samp']], sep='')
fam <- read_fam(snakemake@input[['fam']])
# declare all samples in one group, not using this info
group<-rep(1, nrow(fam))
# get index of related samples
related_samp_which <-which(fam$fam %in% related_samp$V1)
#calculate pca with related subjects projected
sm.pca <- smart_pca(snp_data = geno_loc, pc_axes = 10, sample_group= group, sample_project = related_samp_which, pc_project =1:10)
# format and write out PCs
eigen_vec<- as.data.frame(sm.pca$pca.sample_coordinates)
eigen_vec$ID<- fam$id
eigen_vec <- eigen_vec %>% select(-c(Group,Class))
eigen_vec <- eigen_vec %>% relocate(ID)
write.table(eigen_vec, snakemake@output[['eigenvec']], sep='\t', row.names=FALSE)

