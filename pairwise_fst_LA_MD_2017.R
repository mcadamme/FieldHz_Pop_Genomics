#Script to get overall divergence between MD and LA populations
#07252020 MF


library(adegenet); library(vcfR); library(hierfstat); library(bootstrap)

setwd("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/two_state_2017_samps/")

data <- read.vcfR("Hzea_2017_samps_variantsonly.vcf.recode.vcf")
data_genind <- vcfR2genind(data) #converts file format
data_genind@pop <- as.factor(c(rep("LA", each = 11), rep("MD", each = 16)))

All_hier <- genind2hierfstat(data_genind, pop = data_genind@pop)
sum_data_All_hier <- basic.stats(All_hier)


#getting FST plus bootstrapped confint
print(sum_data_All_hier$overall[7])

Fst_vals <- as.vector(sum_data_All_hier$perloc[,7])

FST_confint <- bcanon(Fst_vals,1000,mean)   
print(FST_confint$confpoints)

