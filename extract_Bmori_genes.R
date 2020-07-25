#Script to extract Bmori genes between outliers on Chr 13
#07232020 MF

library(ape)

gff_Bm <- read.gff("/home/megan/Genomes/Bmori_GCF_000151625.1_ASM15162v1_genomic.gff")
head(gff_Bm)

#getting Chr13 only
Chr13_gff_Bm <- subset(gff_Bm, seqid == "NW_004582008.1" & type == "CDS")

#getting regions of synteny with H. zea outliers
Bm_Chr13_syn <- subset(Chr13_gff_Bm, start > 10900000 & start < 13300000)

write.table(Bm_Chr13_syn, file = "/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/SyntenicToSweepRegion_Bmori.txt")

#getting the list of uniques
unique_genes <- unique(Bm_Chr13_syn$attributes)

write.table(unique_genes, file = "/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/unique_genes.txt")
