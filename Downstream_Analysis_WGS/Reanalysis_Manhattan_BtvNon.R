#Here is the script I used to generate a 10kb windowed manhattan plot for comparison of Bt and non-Bt crops
#04122021 MF

#generating CMplots and getting scaffolds with outliers.

library(CMplot); library(dplyr); library(tidyr)

setwd("~/Desktop/Hz_fieldColl_pop_gen/Reanalysis_PNAS")

#loading zea superscaffolds
Full_Ord_Scafs <- read.table(file = "Hzea_superScaf_genome.txt", header = T)

sig_BAP11A1_DD <- read.csv("sig_BAP11A1_DD.csv", header = T)
sig_BAP11A1_CL <- read.csv("sig_BAP11A1_CL.csv", header = T)
sig_DEO9A1_DD <- read.csv("sig_DEO9A1_DD.csv", header = T)
sig_DEO9A1_CL <- read.csv("sig_DEO9A1_CL.csv", header = T)

#prep for merge
sig_BAP11A1_DD$ScafPos <- paste0(sig_BAP11A1_DD$scaffold,"_", sig_BAP11A1_DD$snp_bin)
sig_BAP11A1_CL$ScafPos <- paste0(sig_BAP11A1_CL$scaffold,"_", sig_BAP11A1_CL$snp_bin)
sig_DEO9A1_DD$ScafPos <- paste0(sig_DEO9A1_DD$scaffold,"_", sig_DEO9A1_DD$snp_bin)
sig_DEO9A1_CL$ScafPos <- paste0(sig_DEO9A1_CL$scaffold,"_", sig_DEO9A1_CL$snp_bin)

#identifying and removing scaffolds shared by both DD and CL for each family
intersect(unique(sig_BAP11A1_CL$scaffold), unique(sig_BAP11A1_DD$scaffold)) #"KZ118207.1" "KZ116810.1"
intersect(unique(sig_DEO9A1_CL$scaffold), unique(sig_DEO9A1_DD$scaffold))

sig_BAP11A1_DD <- subset(sig_BAP11A1_DD, scaffold != "KZ118207.1" & scaffold != "KZ116810.1")

#overlap between families?  Not removing
intersect(unique(sig_BAP11A1_DD$scaffold), unique(sig_DEO9A1_DD$scaffold))#"KZ118135.1" "KZ118249.1"

#choosing sigLev
sig_BAP11A1_DD <- subset(sig_BAP11A1_DD, sigLev == 0.01)
sig_DEO9A1_DD <- subset(sig_DEO9A1_DD, sigLev == 0.01)

##### 2017 Bt and Non-Bt comparison - 10kb #######

wcFST_Bt_unfilt <- read.table("BtandNonBt_Hzea_variantsonly_wcST_10kb_all.smoothed", header = F)
names(wcFST_Bt_unfilt) <- c("Scaf", "WinStart", "WinStop", "NumSnps", "wcFST")
wcFST_Bt_unfilt$Scaf <- as.character(wcFST_Bt_unfilt$Scaf)

wcFST_Bt <- subset(wcFST_Bt_unfilt, NumSnps > 10)

#getting significance threshold based on ztransformation
wcFST_Bt$ztrans <- scale(wcFST_Bt$wcFST, center = TRUE, scale = TRUE) #ztransformation

hist(wcFST_Bt$ztrans, breaks = 100)


#merging ordered scafs with SNPs for plot
forPlot_Bt <- data.frame(cbind(wcFST_Bt$Scaf, wcFST_Bt$WinStart, wcFST_Bt$wcFST, wcFST_Bt$ztrans))
names(forPlot_Bt) <- c("Scaf", "WinStart", "wcFST", "ztrans")
forPlot_Bt$ztrans <- as.numeric(as.character(forPlot_Bt$ztrans)) 

merged_Bt <- merge(Full_Ord_Scafs, forPlot_Bt, by = "Scaf")
merged_Bt$WinStart <- as.numeric(as.character(merged_Bt$WinStart))
merged_Bt <- merged_Bt[order(merged_Bt[,5], merged_Bt[,6], merged_Bt[,3], merged_Bt[,1], merged_Bt[,7]),]


merged_Bt$wcFST <- as.numeric(as.character(merged_Bt$wcFST))
merged_Bt$SnpName <- seq(1:NROW(merged_Bt))
merged_Bt$ScafPos <- paste0(merged_Bt$Scaf,"_",merged_Bt$WinStart)

#getting highlighted SNPs
hi_wcFST0 <- subset(merged_Bt, ztrans > 6) #subsetting out any window with an FST greater than 6 SDs from mean.
print(min(hi_wcFST0$wcFST))
print(mean(wcFST_Bt$wcFST))
print(median(wcFST_Bt$wcFST))

sig_BAP11A1_DD1 <- merge(sig_BAP11A1_DD, merged_Bt, by = "ScafPos")
sig_BAP11A1_DD1$SnpName <- as.numeric(sig_BAP11A1_DD1$SnpName)

sig_DEO9A1_DD1 <- merge(sig_DEO9A1_DD, merged_Bt, by = "ScafPos")
sig_DEO9A1_DD1$SnpName <- as.numeric(sig_DEO9A1_DD1$SnpName)


#highlighting armigera scaffolds with QTLs - cry1
min_HiCry1ByScaf <- as.numeric(as.character(as.vector(tapply(sig_BAP11A1_DD1$SnpName, sig_BAP11A1_DD1$Chr, min))))
max_HiCry1ByScaf <- as.numeric(as.character(as.vector(tapply(sig_BAP11A1_DD1$SnpName, sig_BAP11A1_DD1$Chr, max))))

minMax_cry1 <- cbind(min_HiCry1ByScaf, max_HiCry1ByScaf)
str(minMax_cry1)
names(minMax_cry1) <- c("Min", "Max")

#loop gets the windows
cry1_SnpNames <- vector()
for (row in 1:nrow(minMax_cry1)){
  SNPseq <- (seq(from = minMax_cry1[row,1], to = minMax_cry1[row,2]))
  cry1_SnpNames <- c(cry1_SnpNames, SNPseq)
  
}

#then subset by significant fst value
cry1_SnpNames_DF <- data.frame(cry1_SnpNames)
Cry1FST <- merge(merged_Bt, cry1_SnpNames_DF, by.x = "SnpName", by.y = "cry1_SnpNames")
hiCry1FST <- subset(Cry1FST, ztrans > 6)

write.table(hiCry1FST, "hiCry1FST_BtvNon.txt", row.names = F)


#highlighting all SNPs within QTL region - cry1 + cry2
min_HiCry2ByScaf <- as.numeric(as.character(as.vector(tapply(sig_DEO9A1_DD1$SnpName, sig_DEO9A1_DD1$Chr, min))))
max_HiCry2ByScaf <- as.numeric(as.character(as.vector(tapply(sig_DEO9A1_DD1$SnpName, sig_DEO9A1_DD1$Chr, max))))

minMax_cry2 <- cbind(min_HiCry2ByScaf, max_HiCry2ByScaf)
str(minMax_cry2)
names(minMax_cry2) <- c("Min", "Max")

#loop gets the windows
cry2_SnpNames <- vector()
for (row in 1:nrow(minMax_cry2)){
  SNPseq <- (seq(from = minMax_cry2[row,1], to = minMax_cry2[row,2]))
  cry2_SnpNames <- c(cry2_SnpNames, SNPseq)
  
}

#then subset by significant fst value
cry2_SnpNames_DF <- data.frame(cry2_SnpNames)
Cry2FST <- merge(merged_Bt, cry2_SnpNames_DF, by.x = "SnpName", by.y = "cry2_SnpNames")
hiCry2FST <- subset(Cry2FST, ztrans > 6)

write.table(hiCry2FST, "hiCry2FST_BtvNon.txt", row.names = F)

#no Obs on Chr44 after filter, so needed a work-around for the loop
merged_Bt$Chr_fac <- as.factor(merged_Bt$Chr)
levels(merged_Bt$Chr_fac)[levels(merged_Bt$Chr_fac)=="45"] <- "44"
levels(merged_Bt$Chr_fac)[levels(merged_Bt$Chr_fac)=="46"] <- "45"

#adding artificial continuous window
with_artWin <- data.frame()

for (i in seq_along(levels(merged_Bt$Chr_fac))){
  print(i)
  sub_chrom <- subset(merged_Bt, Chr_fac == i)
  print(nrow(sub_chrom))
  artWin <- seq(0, (nrow(sub_chrom)-1)*10000, by = 10000)
  comb <- data.frame(cbind(sub_chrom, artWin))
  with_artWin <- rbind(with_artWin, comb)
}

with_artWin$wcFST <- as.numeric(as.character(with_artWin$wcFST))
merged_Bt_forPlot <- with_artWin[,c(10,5,13,8)]
names(merged_Bt_forPlot) <- c("SnpName", "Chr", "artWin", "wcFST")

#Manhattan plot function

#sig FST
CMplot(merged_Bt_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST0$wcFST)),
       chr.den.col="pink", file="jpg", memo="BtvNon_10kb_wcFST", dpi=300, threshold = min(hi_wcFST0$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hi_wcFST0$SnpName)

#sig Cry1
CMplot(merged_Bt_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST0$wcFST)),
       chr.den.col="pink", file="jpg", memo="BtvNon_10kb_cry1", dpi=300, threshold = min(hi_wcFST0$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hiCry1FST$SnpName)

#sig Cry2
CMplot(merged_Bt_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST0$wcFST)),
       chr.den.col="pink", file="jpg", memo="BtvNon_10kb_cry2", dpi=300, threshold = min(hi_wcFST0$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hiCry2FST$SnpName)


##### Nearby Genes #####
#pulling out the parts of the gff3 file that are important
gff3 <- read.table("~/Genomes/HzOGS2-15205-fixed_note-added.gff3", sep="\t", stringsAsFactors = F)
head(gff3)

gff3_genes <- subset(gff3, V3 == "gene")
ordered_gff3_genes <- gff3_genes[order(gff3_genes[,1], gff3_genes[,4]),]


scaf_names <- read.table("~/Desktop/Hz_fieldColl_pop_gen/Reanalysis_PNAS/scaffold_and_contig_names.txt", header = T)

all_scafs <- unique(scaf_names$CHROM)


#Bt_v_non-Bt from MD
outliers_MD <- rbind(hiCry1FST, hiCry2FST)
outliers_MD_withIDs <- merge(outliers_MD, scaf_names, by.x = "Scaf", by.y = "CHROM")

uniqScafs_MD <- as.character(unique(outliers_MD_withIDs$Scaf_ID))
uniqScafs_MD

all_genes_outliers_MD <- data.frame()

for (item in 1:length(uniqScafs_MD)){
  sample <- uniqScafs_MD[item]
  print(sample)
  Genes_MD <- ordered_gff3_genes[grepl(sample, ordered_gff3_genes$V1), ]
  all_genes_outliers_MD <- rbind(Genes_MD,all_genes_outliers_MD)
}

outliers_MD_withAnnot <- merge(outliers_MD_withIDs, all_genes_outliers_MD, by.x = "Scaf_ID", by.y = "V1")
str(outliers_MD_withAnnot)

sub_outliers_MD_withAnnot <- subset(outliers_MD_withAnnot, V5 > (WinStart - 50000) & V4 < (WinStart + 50000))
for_FileS4_MD <- sub_outliers_MD_withAnnot[,c(1,2,7,15,16,18,20)]
colnames(for_FileS4_MD) <- c("Scaf_ID", "Scaf", "Chr", "start", "stop", "strand", "gene")

write.table(for_FileS4_MD %>% distinct(), file = "Genes_50kb_outliers_BtvNon.txt") 


