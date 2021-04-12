#Used this script to calculate genomic divergence (average wcFST) across candidate genes for Bt resistance.
#01222020

library(CMplot)

setwd("~/Desktop/Hz_fieldColl_pop_gen/Reanalysis_PNAS/cand_gene")

#loading data files to calculate average for each gene
alp <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/KZ117832.1_only/KZ117832.1_wcFst_2002_2017", header = F)
apn1 <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/KZ118301.1_only/KZ118301.1_wcFst_2002_2017", header = F)
apn4 <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/KZ118301.1_only/KZ118301.1_wcFst_2002_2017", header = F)
abcA2 <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/KZ118207.1_only/KZ118207.1_wcFst_2002_2017", header = F)
abcC2 <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/KZ118297.1_only/KZ118297.1_wcFst_2002_2017", header = F)
abcG1 <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/KZ118424.1_only/KZ118424.1_wcFst_2002_2017", header = F)
calp <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/KZ117563.1_only/KZ117563.1_wcFst_2002_2017", header = F)
cad2 <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/KZ118195.1_only/KZ118195.1_wcFst_2002_2017", header = F)
cad86c <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1//KZ117463.1_only/KZ117463.1_wcFst_2002_2017", header = F)
map4K4 <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/KZ118817.1_only/KZ118817.1_wcFst_2002_2017", header = F)
tspan1 <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/KZ118424.1_only/KZ118424.1_wcFst_2002_2017", header = F)

carbQ <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/KZ118395.1_only/KZ118395.1_wcFst_2002_2017", header = F)
cyp333b <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/KZ117131.1_only/KZ117131.1_wcFst_2002_2017", header = F)
venpep <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/KZ117237.1_only/KZ117237.1_wcFst_2002_2017", header = F)

#getting gene regions
alp_gene <- subset(alp, V2 > 1946 & V2 < 16317)
mean(alp_gene$V5)
nrow(alp_gene)

apn1_gene <- subset(apn1, V2 > 154624 & V2 < 162221)
mean(apn1_gene$V5)
nrow(apn1_gene)

apn4_gene <- subset(apn4, V2 > 167611 & V2 < 174331)
mean(apn4_gene$V5)
nrow(apn4_gene)

abcA2_gene <- subset(abcA2, V2 > 174568 & V2 < 191566)
mean(abcA2_gene$V5)
nrow(abcA2_gene)

abcC2_gene <- subset(abcC2, V2 > 98882 & V2 < 114139)
mean(abcC2_gene$V5)
nrow(abcC2_gene)

abcG1_gene <- subset(abcG1, V2 > 8580 & V2 < 37181)
mean(abcG1_gene$V5)
nrow(abcG1_gene)

calp_gene <- subset(calp, V2 > 270372 & V2 < 288777)
mean(calp_gene$V5)
nrow(calp_gene)

cad2_gene <- subset(cad2, V2 > 87158 & V2 < 103404)
mean(cad2_gene$V5)
nrow(cad2_gene)

cad86c_gene <- subset(cad86c, V2 > 550935 & V2 < 603715)
mean(cad86c_gene$V5)
nrow(cad86c_gene)

map4K4_gene <- subset(map4K4, V2 > 92551 & V2 < 104845)
mean(map4K4_gene$V5)
nrow(map4K4_gene)

tspan1_gene <- subset(tspan1, V2 > 105935 & V2 < 117025)
mean(tspan1_gene$V5)
nrow(tspan1_gene)

#from blast
carbQ_gene1 <- subset(carbQ, V2 > 8573 & V2 < 10730)
mean(carbQ_gene1$V5)
nrow(carbQ_gene1)

#from annotation
carbQ_gene2 <- subset(carbQ, V2 > 151736 & V2 < 160706)
mean(carbQ_gene2$V5)
nrow(carbQ_gene2)

cyp333b_gene <- subset(cyp333b, V2 > 32819 & V2 < 39355)
mean(cyp333b_gene$V5)
nrow(cyp333b_gene)

venpep_gene <- subset(venpep, V2 > 266124 & V2 < 298863)
mean(venpep_gene$V5)
nrow(venpep_gene)

#loading data to plot sliding window averaged fst - 10kb windows with 1kb step size

alp_win <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/KZ117832.1_only/KZ117832.1_wcFst_2002_2017.smoothed", header = F)
apn1_win <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/KZ118301.1_only/KZ118301.1_wcFst_2002_2017.smoothed", header = F)
apn4_win <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/KZ118301.1_only/KZ118301.1_wcFst_2002_2017.smoothed", header = F)
abcA2_win <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/KZ118207.1_only/KZ118207.1_wcFst_2002_2017.smoothed", header = F)
abcC2_win <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/KZ118297.1_only/KZ118297.1_wcFst_2002_2017.smoothed", header = F)
abcG1_win <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/KZ118424.1_only/KZ118424.1_wcFst_2002_2017.smoothed", header = F)
calp_win <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/KZ117563.1_only/KZ117563.1_wcFst_2002_2017.smoothed", header = F)
cad2_win <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/KZ118195.1_only/KZ118195.1_wcFST_2002_2017.smoothed", header = F)
cad86c_win <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/KZ117463.1_only/KZ117463.1_wcFst_2002_2017.smoothed", header = F)
map4K4_win <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/KZ118817.1_only/KZ118817.1_wcFst_2002_2017.smoothed", header = F)
tspan1_win <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/KZ118424.1_only/KZ118424.1_wcFst_2002_2017.smoothed", header = F)

carbQ_win <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/KZ118395.1_only/KZ118395.1_wcFst_2002_2017.smoothed", header = F)
cyp333b_win <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/KZ117131.1_only/KZ117131.1_wcFst_2002_2017.smoothed", header = F)
venpep_win <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/KZ117237.1_only/KZ117237.1_wcFst_2002_2017.smoothed", header = F)



#Manhattan FST plot function - for candidate genes
window_name <- function(dataset_win, dataset_gene){
  win <- seq(1:(nrow(dataset_win)))
  print(win)
  new_data <- cbind(win,dataset_win)
  new_data <- data.frame(new_data)
  contig_name <- as.character(dataset_win[1,1])
  str(new_data)
  
  gene <- subset(new_data, V2 > min(dataset_gene$V2) & V2 < max(dataset_gene$V2))
  head(gene)
  HiLite <- as.numeric(as.character(gene$win))
  print(HiLite)
  
  for_plot <- data.frame(cbind(win, new_data$V1, new_data$V2, new_data$V5))
  print(for_plot)
  
  CMplot(for_plot, plot.type="m", type = "l", r=1.6, cir.legend=TRUE, col = c("grey30"), cex = 0.8, ylim = c(-0.05,0.25),
         file="jpg", memo="", dpi=300, threshold = NULL, LOG10 = F, highlight = HiLite, ylab = "FST", xlab = contig_name)
  
  
}

#dataset names get changed here - then must go rename file
window_name(tspan1_win, tspan1_gene)


#Manhattan FST - for combined carboxyQ plot given that there appear to be two genes on this scaffold.
win <- seq(1:(nrow(carbQ_win)))
print(win)

new_data <- cbind(win,carbQ_win)
new_data <- data.frame(new_data)
contig_name <- as.character(carbQ_win[1,1])
str(new_data)

gene1 <- subset(new_data, V2 > min(carbQ_gene1$V2) & V2 < max(carbQ_gene1$V2))
head(gene1)

gene2 <- subset(new_data, V2 > min(carbQ_gene2$V2) & V2 < max(carbQ_gene2$V2))
head(gene2)

both_genes <- rbind(gene1,gene2) 
HiLite <- as.numeric(as.character(both_genes$win))
print(HiLite)

for_plot <- data.frame(cbind(win, new_data$V1, new_data$V2, new_data$V5))
print(for_plot)

CMplot(for_plot, plot.type="m", type = "l", r=1.6, cir.legend=TRUE, col = c("grey30"), cex = 0.8, ylim = c(-0.05,0.25),
       file="jpg", memo="", dpi=300, threshold = NULL, LOG10 = F, highlight = HiLite, ylab = "FST", xlab = "Scaffold 569")



#Looking at cpq, cyp scaffolds for Bt and NonBt samples
carbQ_Bt <- read.table("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/BtandNonBt_alignmentFiles/WGRS_mpileupANDvcftools_output/KZ118395.1_only/BtandNonBt_Hzea_variantsonly_wcST_10kb_KZ118395.1.smoothed", header = F)
#cyp333b <- read.table("./KZ117131.1_only/KZ117131.1_wcFst_2002_2017", header = F)

win <- seq(1:(nrow(carbQ_Bt)))
print(win)

new_data <- cbind(win,carbQ_Bt)
new_data <- data.frame(new_data)
contig_name <- as.character(carbQ_Bt[1,1])
str(new_data)

gene1 <- subset(new_data, V2 > min(carbQ_gene1$V2) & V2 < max(carbQ_gene1$V2))
head(gene1)

gene2 <- subset(new_data, V2 > min(carbQ_gene2$V2) & V2 < max(carbQ_gene2$V2))
head(gene2)

both_genes <- rbind(gene1,gene2) 
HiLite <- as.numeric(as.character(both_genes$win))
print(HiLite)

for_plot <- data.frame(cbind(win, new_data$V1, new_data$V2, new_data$V5))
print(for_plot)

CMplot(for_plot, plot.type="m", type = "l", r=1.6, cir.legend=TRUE, col = c("grey30"), cex = 0.8, ylim = c(-0.05,0.25),
       file="jpg", memo="BtandNonBtScaf569", dpi=300, threshold = NULL, LOG10 = F, highlight = HiLite, ylab = "FST", xlab = "Scaffold 569")
