#Used this script to calculate genomic divergence (average wcFST) across candidate genes for Bt resistance.
#01222020

library(CMplot)

setwd("/media/megan/New Volume/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/")

#loading data files to calculate average for each gene
alp <- read.table("./KZ117832.1_only/KZ117832.1_wcFst_2002_2017", header = F)
apn1 <- read.table("./KZ118301.1_only/KZ118301.1_wcFst_2002_2017", header = F)
apn4 <- read.table("./KZ118301.1_only/KZ118301.1_wcFst_2002_2017", header = F)
abcA2 <- read.table("./KZ118207.1_only/KZ118207.1_wcFst_2002_2017", header = F)
abcC2 <- read.table("./KZ118297.1_only/KZ118297.1_wcFst_2002_2017", header = F)
abcG1 <- read.table("./KZ118424.1_only/KZ118424.1_wcFst_2002_2017", header = F)
calp <- read.table("./KZ117563.1_only/KZ117563.1_wcFst_2002_2017", header = F)
cad2 <- read.table("./KZ118195.1_only/KZ118195.1_wcFst_2002_2017", header = F)
cad86c <- read.table("./KZ117463.1_only/KZ117463.1_wcFst_2002_2017", header = F)
map4K4 <- read.table("./KZ118817.1_only/KZ118817.1_wcFst_2002_2017", header = F)
tspan1 <- read.table("./KZ118424.1_only/KZ118424.1_wcFst_2002_2017", header = F)

carbQ <- read.table("./KZ118395.1_only/KZ118395.1_wcFst_2002_2017", header = F)
cyp333b <- read.table("./KZ117131.1_only/KZ117131.1_wcFst_2002_2017", header = F)
venpep <- read.table("./KZ117237.1_only/KZ117237.1_wcFst_2002_2017", header = F)

#getting gene regions
alp_gene <- subset(alp, V2 > 1946 & V2 < 16347)
mean(alp_gene$V5)
nrow(alp_gene)

apn1_gene <- subset(apn1, V2 > 140429 & V2 < 161526)
mean(apn1_gene$V5)
nrow(apn1_gene)

apn4_gene <- subset(apn4, V2 > 168017 & V2 < 174495)
mean(apn4_gene$V5)
nrow(apn4_gene)

abcA2_gene <- subset(abcA2, V2 > 171103 & V2 < 192554)
mean(abcA2_gene$V5)
nrow(abcA2_gene)

abcC2_gene <- subset(abcC2, V2 > 98882 & V2 < 114139)
mean(abcC2_gene$V5)
nrow(abcC2_gene)

abcG1_gene <- subset(abcG1, V2 > 2680 & V2 < 24508)
mean(abcG1_gene$V5)
nrow(abcG1_gene)

calp_gene <- subset(calp, V2 > 271234 & V2 < 288320)
mean(calp_gene$V5)
nrow(calp_gene)

cad2_gene <- subset(cad2, V2 > 87164 & V2 < 95177)
mean(cad2_gene$V5)
nrow(cad2_gene)

cad86c_gene <- subset(cad86c, V2 > 550935 & V2 < 603715)
mean(cad86c_gene$V5)
nrow(cad86c_gene)

map4K4_gene <- subset(map4K4, V2 > 92551 & V2 < 104845)
mean(map4K4_gene$V5)
nrow(map4K4_gene)

tspan1_gene <- subset(tspan1, V2 > 106132 & V2 < 117025)
mean(tspan1_gene$V5)
nrow(tspan1_gene)

carbQ_gene <- subset(carbQ, V2 > 151736 & V2 < 160706)
mean(carbQ_gene$V5)
nrow(carbQ_gene)

cyp333b_gene <- subset(cyp333b, V2 > 32819 & V2 < 39355)
mean(cyp333b_gene$V5)
nrow(cyp333b_gene)

venpep_gene <- subset(venpep, V2 > 266124 & V2 < 298863)
mean(venpep_gene$V5)
nrow(venpep_gene)

#loading data to plot sliding window averaged fst - 10kb windows with 1kb step size

alp_win <- read.table("./KZ117832.1_only/KZ117832.1_wcFst_2002_2017.smoothed", header = F)
apn1_win <- read.table("./KZ118301.1_only/KZ118301.1_wcFst_2002_2017.smoothed", header = F)
apn4_win <- read.table("./KZ118301.1_only/KZ118301.1_wcFst_2002_2017.smoothed", header = F)
abcA2_win <- read.table("./KZ118207.1_only/KZ118207.1_wcFst_2002_2017.smoothed", header = F)
abcC2_win <- read.table("./KZ118297.1_only/KZ118297.1_wcFst_2002_2017.smoothed", header = F)
abcG1_win <- read.table("./KZ118424.1_only/KZ118424.1_wcFst_2002_2017.smoothed", header = F)
calp_win <- read.table("./KZ117563.1_only/KZ117563.1_wcFst_2002_2017.smoothed", header = F)
cad2_win <- read.table("./KZ118195.1_only/KZ118195.1_wcFST_2002_2017.smoothed", header = F)
cad86c_win <- read.table("./KZ117463.1_only/KZ117463.1_wcFst_2002_2017.smoothed", header = F)
map4K4_win <- read.table("./KZ118817.1_only/KZ118817.1_wcFst_2002_2017.smoothed", header = F)
tspan1_win <- read.table("./KZ118424.1_only/KZ118424.1_wcFst_2002_2017.smoothed", header = F)

carbQ_win <- read.table("./KZ118395.1_only/KZ118395.1_wcFst_2002_2017.smoothed", header = F)
cyp333b_win <- read.table("./KZ117131.1_only/KZ117131.1_wcFst_2002_2017.smoothed", header = F)
venpep_win <- read.table("./KZ117237.1_only/KZ117237.1_wcFst_2002_2017.smoothed", header = F)



#Manhattan FST plot function
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
   

