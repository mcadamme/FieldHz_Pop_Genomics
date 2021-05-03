#Here is the script I used to generate new Manhattan plots with 10kb window sizes and add QTL
#03162021 MF

#generating CMplots and getting scaffolds with outliers.

library(CMplot); library(dplyr); library(tidyr)

setwd("~/Desktop/Hz_fieldColl_pop_gen/Reanalysis_PNAS")

#loading zea superscaffolds
Full_Ord_Scafs <- read.table(file = "Hzea_superScaf_genome.txt", header = T)

#loading QTL SNPs
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


##### 2002-2017 comparison - 10kb #######

wcFST_2002_2017_unfilt <- read.table("2002and2017_10kb_wcFST_all.smoothed", header = F)
names(wcFST_2002_2017_unfilt) <- c("Scaf", "WinStart", "WinStop", "NumSnps", "wcFST")
wcFST_2002_2017_unfilt$Scaf <- as.character(wcFST_2002_2017_unfilt$Scaf)

wcFST_2002_2017 <- subset(wcFST_2002_2017_unfilt, NumSnps > 10)

#getting significance threshold based on ztransformation
wcFST_2002_2017$ztrans <- scale(wcFST_2002_2017$wcFST, center = TRUE, scale = TRUE) #ztransformation

hist(wcFST_2002_2017$ztrans, breaks = 100)


#merging ordered scafs with SNPs for plot
forPlot_2002_2017 <- data.frame(cbind(wcFST_2002_2017$Scaf, wcFST_2002_2017$WinStart, wcFST_2002_2017$wcFST, wcFST_2002_2017$ztrans))
names(forPlot_2002_2017) <- c("Scaf", "WinStart", "wcFST", "ztrans")
forPlot_2002_2017$ztrans <- as.numeric(as.character(forPlot_2002_2017$ztrans)) 

merged_2002_2017 <- merge(Full_Ord_Scafs, forPlot_2002_2017, by = "Scaf")
merged_2002_2017$WinStart <- as.numeric(as.character(merged_2002_2017$WinStart))
merged_2002_2017 <- merged_2002_2017[order(merged_2002_2017 [,5], merged_2002_2017 [,6], merged_2002_2017 [,3], 
                                     merged_2002_2017 [,1], merged_2002_2017 [,7]),]


merged_2002_2017$wcFST <- as.numeric(as.character(merged_2002_2017$wcFST))
merged_2002_2017$SnpName <- seq(1:NROW(merged_2002_2017))
merged_2002_2017$ScafPos <- paste0(merged_2002_2017$Scaf,"_",merged_2002_2017$WinStart)

#getting highlighted SNPs
hi_wcFST0 <- subset(merged_2002_2017, ztrans > 6) #subsetting out any window with an FST greater than 6 SDs from mean.
nrow(hi_wcFST0)#gives number of sig SNPs
length(unique(hi_wcFST0$Scaf)) #gives number of scaffolds
length(unique(hi_wcFST0$Chr)) #gives number of Chromosomes

print(min(hi_wcFST0$wcFST))
print(max(hi_wcFST0$wcFST))
print(mean(wcFST_2002_2017$wcFST))
print(median(wcFST_2002_2017$wcFST))

sig_BAP11A1_DD1 <- merge(sig_BAP11A1_DD, merged_2002_2017, by = "ScafPos")
sig_BAP11A1_DD1$SnpName <- as.numeric(sig_BAP11A1_DD1$SnpName)

sig_DEO9A1_DD1 <- merge(sig_DEO9A1_DD, merged_2002_2017, by = "ScafPos")
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
Cry1FST <- merge(merged_2002_2017, cry1_SnpNames_DF, by.x = "SnpName", by.y = "cry1_SnpNames")
hiCry1FST <- subset(Cry1FST, ztrans > 6)

write.table(hiCry1FST, "hiCry1FST2002to2017.txt", row.names = F)


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
Cry2FST <- merge(merged_2002_2017, cry2_SnpNames_DF, by.x = "SnpName", by.y = "cry2_SnpNames")
hiCry2FST <- subset(Cry2FST, ztrans > 6)

write.table(hiCry2FST, "hiCry2FST2002to2017.txt", row.names = F)

#no Obs on Chr44 after filter, so needed a work-around for the loop
merged_2002_2017$Chr_fac <- as.factor(merged_2002_2017$Chr)
levels(merged_2002_2017$Chr_fac)[levels(merged_2002_2017$Chr_fac)=="45"] <- "44"
levels(merged_2002_2017$Chr_fac)[levels(merged_2002_2017$Chr_fac)=="46"] <- "45"

#adding artificial continuous window
with_artWin <- data.frame()

for (i in seq_along(levels(merged_2002_2017$Chr_fac))){
  print(i)
  sub_chrom <- subset(merged_2002_2017, Chr_fac == i)
  print(nrow(sub_chrom))
  artWin <- seq(0, (nrow(sub_chrom)-1)*10000, by = 10000)
  comb <- data.frame(cbind(sub_chrom, artWin))
  with_artWin <- rbind(with_artWin, comb)
}

with_artWin$wcFST <- as.numeric(as.character(with_artWin$wcFST))
merged_2002_2017_forPlot <- with_artWin[,c(10,5,13,8)]
names(merged_2002_2017_forPlot) <- c("SnpName", "Chr", "artWin", "wcFST")

#2002-2017 Manhattan plot function

#sig FST
CMplot(merged_2002_2017_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST0$wcFST)),
       chr.den.col="pink", file="jpg", memo="2002and2017_10kb_wcFST", dpi=300, threshold = min(hi_wcFST0$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hi_wcFST0$SnpName)

#sig Cry1
CMplot(merged_2002_2017_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST0$wcFST)),
       chr.den.col="pink", file="jpg", memo="2002and2017_10kb_cry1", dpi=300, threshold = min(hi_wcFST0$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hiCry1FST$SnpName)

#sig Cry2
CMplot(merged_2002_2017_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST0$wcFST)),
       chr.den.col="pink", file="jpg", memo="2002and2017_10kb_cry2", dpi=300, threshold = min(hi_wcFST0$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hiCry2FST$SnpName)


##### 2002-2012 comparison - 10kb #####
wcFST_2002_2012_unfilt <- read.table("2002and2012_10kb_wcFST_all.smoothed", header = F)
names(wcFST_2002_2012_unfilt) <- c("Scaf", "WinStart", "WinStop", "NumSnps", "wcFST")
wcFST_2002_2012_unfilt$Scaf <- as.character(wcFST_2002_2012_unfilt$Scaf)

wcFST_2002_2012 <- subset(wcFST_2002_2012_unfilt, NumSnps > 10)

#getting significance threshold based on ztransformation
wcFST_2002_2012$ztrans <- scale(wcFST_2002_2012$wcFST, center = TRUE, scale = TRUE) #ztransformation

hist(wcFST_2002_2012$ztrans, breaks = 100)


#merging ordered scafs with SNPs for plot
forPlot_2002_2012 <- data.frame(cbind(wcFST_2002_2012$Scaf, wcFST_2002_2012$WinStart, wcFST_2002_2012$wcFST, wcFST_2002_2012$ztrans))
names(forPlot_2002_2012) <- c("Scaf", "WinStart", "wcFST", "ztrans")
forPlot_2002_2012$ztrans <- as.numeric(as.character(forPlot_2002_2012$ztrans)) 

merged_2002_2012 <- merge(Full_Ord_Scafs, forPlot_2002_2012, by = "Scaf")
merged_2002_2012$WinStart <- as.numeric(as.character(merged_2002_2012$WinStart))
merged_2002_2012 <- merged_2002_2012[order(merged_2002_2012 [,5], merged_2002_2012 [,6], merged_2002_2012 [,3], 
                                           merged_2002_2012 [,1], merged_2002_2012 [,7]),]


merged_2002_2012$wcFST <- as.numeric(as.character(merged_2002_2012$wcFST))
merged_2002_2012$SnpName <- seq(1:NROW(merged_2002_2012))
merged_2002_2012$ScafPos <- paste0(merged_2002_2012$Scaf,"_",merged_2002_2012$WinStart)

#getting highlighted SNPs
hi_wcFST1 <- subset(merged_2002_2012, ztrans > 6) #subsetting out any window with an FST greater than 6 SDs from mean.
nrow(hi_wcFST1)#gives number of sig SNPs
length(unique(hi_wcFST1$Scaf)) #gives number of scaffolds
length(unique(hi_wcFST1$Chr)) #gives number of Chromosomes

print(min(hi_wcFST1$wcFST))
print(mean(wcFST_2002_2012$wcFST))
print(median(wcFST_2002_2012$wcFST))

sig_BAP11A1_DD2 <- merge(sig_BAP11A1_DD, merged_2002_2012, by = "ScafPos")
sig_BAP11A1_DD2$SnpName <- as.numeric(sig_BAP11A1_DD2$SnpName)

sig_DEO9A1_DD2 <- merge(sig_DEO9A1_DD, merged_2002_2012, by = "ScafPos")
sig_DEO9A1_DD2$SnpName <- as.numeric(sig_DEO9A1_DD2$SnpName)


#highlighting armigera scaffolds with QTLs - cry1
min_HiCry1ByScaf <- as.numeric(as.character(as.vector(tapply(sig_BAP11A1_DD2$SnpName, sig_BAP11A1_DD2$Chr, min))))
max_HiCry1ByScaf <- as.numeric(as.character(as.vector(tapply(sig_BAP11A1_DD2$SnpName, sig_BAP11A1_DD2$Chr, max))))

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
Cry1FST <- merge(merged_2002_2012, cry1_SnpNames_DF, by.x = "SnpName", by.y = "cry1_SnpNames")
hiCry1FST1 <- subset(Cry1FST, ztrans > 6)

write.table(hiCry1FST1, "hiCry1FST12002to2012.txt", row.names = F)

#highlighting all SNPs within QTL region - cry1 + cry2
min_HiCry2ByScaf <- as.numeric(as.character(as.vector(tapply(sig_DEO9A1_DD2$SnpName, sig_DEO9A1_DD2$Chr, min))))
max_HiCry2ByScaf <- as.numeric(as.character(as.vector(tapply(sig_DEO9A1_DD2$SnpName, sig_DEO9A1_DD2$Chr, max))))

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
Cry2FST <- merge(merged_2002_2012, cry2_SnpNames_DF, by.x = "SnpName", by.y = "cry2_SnpNames")
hiCry2FST1 <- subset(Cry2FST, ztrans > 6)

write.table(hiCry2FST1, "hiCry2FST12002to2012.txt", row.names = F)

#no Obs on Chr44 after filter, so needed a work-around for the loop
merged_2002_2012$Chr_fac <- as.factor(merged_2002_2012$Chr)
levels(merged_2002_2012$Chr_fac)[levels(merged_2002_2012$Chr_fac)=="45"] <- "44"
levels(merged_2002_2012$Chr_fac)[levels(merged_2002_2012$Chr_fac)=="46"] <- "45"

#adding artificial continuous window
with_artWin <- data.frame()

for (i in seq_along(levels(merged_2002_2012$Chr_fac))){
  print(i)
  sub_chrom <- subset(merged_2002_2012, Chr_fac == i)
  print(nrow(sub_chrom))
  artWin <- seq(0, (nrow(sub_chrom)-1)*10000, by = 10000)
  comb <- data.frame(cbind(sub_chrom, artWin))
  with_artWin <- rbind(with_artWin, comb)
}

with_artWin$wcFST <- as.numeric(as.character(with_artWin$wcFST))
merged_2002_2012_forPlot <- with_artWin[,c(10,5,13,8)]
names(merged_2002_2012_forPlot) <- c("SnpName", "Chr", "artWin", "wcFST")

#2002-2012 Manhattan plot function

#sig FST
CMplot(merged_2002_2012_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST1$wcFST)),
       chr.den.col="pink", file="jpg", memo="2002and2012_10kb_wcFST", dpi=300, threshold = min(hi_wcFST1$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hi_wcFST1$SnpName)

#sig Cry1
CMplot(merged_2002_2012_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST1$wcFST)),
       chr.den.col="pink", file="jpg", memo="2002and2012_10kb_cry1", dpi=300, threshold = min(hi_wcFST1$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hiCry1FST1$SnpName)

#sig Cry2
CMplot(merged_2002_2012_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST1$wcFST)),
       chr.den.col="pink", file="jpg", memo="2002and2012_10kb_cry2", dpi=300, threshold = min(hi_wcFST1$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hiCry2FST1$SnpName)


##### 2012-2017 comparison - 10kb #####
wcFST_2012_2017_unfilt <- read.table("2012and2017_10kb_wcFST_all.smoothed", header = F)
names(wcFST_2012_2017_unfilt) <- c("Scaf", "WinStart", "WinStop", "NumSnps", "wcFST")
wcFST_2012_2017_unfilt$Scaf <- as.character(wcFST_2012_2017_unfilt$Scaf)

wcFST_2012_2017 <- subset(wcFST_2012_2017_unfilt, NumSnps > 10)

#getting significance threshold based on ztransformation
wcFST_2012_2017$ztrans <- scale(wcFST_2012_2017$wcFST, center = TRUE, scale = TRUE) #ztransformation

hist(wcFST_2012_2017$ztrans, breaks = 100)


#merging ordered scafs with SNPs for plot
forPlot_2012_2017 <- data.frame(cbind(wcFST_2012_2017$Scaf, wcFST_2012_2017$WinStart, wcFST_2012_2017$wcFST, wcFST_2012_2017$ztrans))
names(forPlot_2012_2017) <- c("Scaf", "WinStart", "wcFST", "ztrans")
forPlot_2012_2017$ztrans <- as.numeric(as.character(forPlot_2012_2017$ztrans)) 

merged_2012_2017 <- merge(Full_Ord_Scafs, forPlot_2012_2017, by = "Scaf")
merged_2012_2017$WinStart <- as.numeric(as.character(merged_2012_2017$WinStart))
merged_2012_2017 <- merged_2012_2017[order(merged_2012_2017 [,5], merged_2012_2017 [,6], merged_2012_2017 [,3], 
                                           merged_2012_2017 [,1], merged_2012_2017 [,7]),]


merged_2012_2017$wcFST <- as.numeric(as.character(merged_2012_2017$wcFST))
merged_2012_2017$SnpName <- seq(1:NROW(merged_2012_2017))
merged_2012_2017$ScafPos <- paste0(merged_2012_2017$Scaf,"_",merged_2012_2017$WinStart)

#getting highlighted SNPs
hi_wcFST2 <- subset(merged_2012_2017, ztrans > 6) #subsetting out any window with an FST greater than 6 SDs from mean.
nrow(hi_wcFST2)#gives number of sig SNPs
length(unique(hi_wcFST2$Scaf)) #gives number of scaffolds
length(unique(hi_wcFST2$Chr)) #gives number of Chromosomes

print(min(hi_wcFST2$wcFST))
print(mean(wcFST_2012_2017$wcFST))
print(median(wcFST_2012_2017_unfilt$wcFST))

sig_BAP11A1_DD3 <- merge(sig_BAP11A1_DD, merged_2012_2017, by = "ScafPos")
sig_BAP11A1_DD3$SnpName <- as.numeric(sig_BAP11A1_DD3$SnpName)

sig_DEO9A1_DD3 <- merge(sig_DEO9A1_DD, merged_2012_2017, by = "ScafPos")
sig_DEO9A1_DD3$SnpName <- as.numeric(sig_DEO9A1_DD3$SnpName)

#highlighting armigera scaffolds with QTLs - cry1
min_HiCry1ByScaf <- as.numeric(as.character(as.vector(tapply(sig_BAP11A1_DD3$SnpName, sig_BAP11A1_DD3$Chr, min))))
max_HiCry1ByScaf <- as.numeric(as.character(as.vector(tapply(sig_BAP11A1_DD3$SnpName, sig_BAP11A1_DD3$Chr, max))))

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
Cry1FST <- merge(merged_2012_2017, cry1_SnpNames_DF, by.x = "SnpName", by.y = "cry1_SnpNames")
hiCry1FST2 <- subset(Cry1FST, ztrans > 6)

write.table(hiCry1FST2, "hiCry1FST22012to2017.txt", row.names = F)

#highlighting all SNPs within QTL region - cry1 + cry2
min_HiCry2ByScaf <- as.numeric(as.character(as.vector(tapply(sig_DEO9A1_DD3$SnpName, sig_DEO9A1_DD3$Chr, min))))
max_HiCry2ByScaf <- as.numeric(as.character(as.vector(tapply(sig_DEO9A1_DD3$SnpName, sig_DEO9A1_DD3$Chr, max))))

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
Cry2FST <- merge(merged_2012_2017, cry2_SnpNames_DF, by.x = "SnpName", by.y = "cry2_SnpNames")
hiCry2FST2 <- subset(Cry2FST, ztrans > 6)

write.table(hiCry2FST2, "hiCry2FST22012to2017.txt", row.names = F)

#no Obs on Chr44 after filter, so needed a work-around for the loop
merged_2012_2017$Chr_fac <- as.factor(merged_2012_2017$Chr)
levels(merged_2012_2017$Chr_fac)[levels(merged_2012_2017$Chr_fac)=="45"] <- "44"
levels(merged_2012_2017$Chr_fac)[levels(merged_2012_2017$Chr_fac)=="46"] <- "45"

#adding artificial continuous window
with_artWin <- data.frame()

for (i in seq_along(levels(merged_2012_2017$Chr_fac))){
  print(i)
  sub_chrom <- subset(merged_2012_2017, Chr_fac == i)
  print(nrow(sub_chrom))
  artWin <- seq(0, (nrow(sub_chrom)-1)*10000, by = 10000)
  comb <- data.frame(cbind(sub_chrom, artWin))
  with_artWin <- rbind(with_artWin, comb)
}

with_artWin$wcFST <- as.numeric(as.character(with_artWin$wcFST))

merged_2012_2017_forPlot <- with_artWin[,c(10,5,13,8)]
names(merged_2012_2017_forPlot) <- c("SnpName", "Chr", "artWin", "wcFST")

#2012-2017 Manhattan plot function

#sig FST
CMplot(merged_2012_2017_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST2$wcFST)),
       chr.den.col="pink", file="jpg", memo="2012and2017_10kb_wcFST", dpi=300, threshold = min(hi_wcFST2$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hi_wcFST2$SnpName)

#sig Cry1
CMplot(merged_2012_2017_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST2$wcFST)),
       chr.den.col="pink", file="jpg", memo="2012and2017_10kb_cry1", dpi=300, threshold = min(hi_wcFST2$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hiCry1FST2$SnpName)

#sig Cry2
CMplot(merged_2012_2017_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST2$wcFST)),
       chr.den.col="pink", file="jpg", memo="2012and2017_10kb_cry2", dpi=300, threshold = min(hi_wcFST2$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hiCry2FST2$SnpName)

#What are the scaffolds not represented in this final plot?

scaf_names <- read.table("~/Desktop/Hz_fieldColl_pop_gen/Reanalysis_PNAS/scaffold_and_contig_names.txt", header = T)

all_scafs <- unique(scaf_names$CHROM)
plotted_scafs <- unique(Full_Ord_Scafs$Scaf)
diff <- setdiff(all_scafs, plotted_scafs)
diff_df <- data.frame(diff)

#no candidate genes from Table 1 on these scaffolds - checked by hand.

intersect(diff_df$diff, unique(hi_wcFST0$Scaf))#as of 3/22 this should be 0.
#added "KZ116441.1" "KZ116446.1" "KZ117435.1" to ord_Chroms_final on 3/22/21, as chromsomes 43-45.  
#These could not be mapped to chromosomes, even by hand.  No known genes on them. Very short.


##### Getting genes near outliers #####

#Now pulling out the parts of the gff3 file that are important
gff3 <- read.table("~/Genomes/HzOGS2-15205-fixed_note-added.gff3", sep="\t", stringsAsFactors = F)
head(gff3)

gff3_genes <- subset(gff3, V3 == "gene")
ordered_gff3_genes <- gff3_genes[order(gff3_genes[,1], gff3_genes[,4]),]


#comparisons for pairs of years - all 10kb windows
#2002-2017
outliers_2002and2017 <- rbind(hiCry1FST, hiCry2FST)
outliers_2002and2017_withIDs <- merge(outliers_2002and2017, scaf_names, by.x = "Scaf", by.y = "CHROM")

uniqScafs_2002and2017 <- as.character(unique(outliers_2002and2017_withIDs$Scaf_ID))
uniqScafs_2002and2017

all_genes_outliers_2002and2017 <- data.frame()

for (item in 1:length(uniqScafs_2002and2017)){
  sample <- uniqScafs_2002and2017[item]
  print(sample)
  Genes_2002and2017 <- ordered_gff3_genes[grepl(sample, ordered_gff3_genes$V1), ]
  all_genes_outliers_2002and2017 <- rbind(Genes_2002and2017,all_genes_outliers_2002and2017)
}

outliers_2002and2017_withAnnot <- merge(outliers_2002and2017_withIDs, all_genes_outliers_2002and2017, by.x = "Scaf_ID", by.y = "V1")
str(outliers_2002and2017_withAnnot)

sub_outliers_2002and2017_withAnnot <- subset(outliers_2002and2017_withAnnot, V5 > (WinStart - 50000) & V4 < (WinStart + 50000))
for_FileS4_2002and2017 <- sub_outliers_2002and2017_withAnnot[,c(1,2,7,15,16,18,20)]
colnames(for_FileS4_2002and2017) <- c("Scaf_ID", "Scaf", "Chr", "start", "stop", "strand", "gene")

write.table(for_FileS4_2002and2017 %>% distinct(), file = "Genes_50kb_outliers_2002&2017.txt") 


#2002-2012
outliers_2002and2012 <- rbind(hiCry1FST1, hiCry2FST1)
outliers_2002and2012_withIDs <- merge(outliers_2002and2012, scaf_names, by.x = "Scaf", by.y = "CHROM")

uniqScafs_2002and2012 <- as.character(unique(outliers_2002and2012_withIDs$Scaf_ID))
uniqScafs_2002and2012

all_genes_outliers_2002and2012 <- data.frame()

for (item in 1:length(uniqScafs_2002and2012)){
  sample <- uniqScafs_2002and2012[item]
  print(sample)
  Genes_2002and2012 <- ordered_gff3_genes[grepl(sample, ordered_gff3_genes$V1), ]
  all_genes_outliers_2002and2012 <- rbind(Genes_2002and2012,all_genes_outliers_2002and2012)
}

outliers_2002and2012_withAnnot <- merge(outliers_2002and2012_withIDs, all_genes_outliers_2002and2012, by.x = "Scaf_ID", by.y = "V1")
str(outliers_2002and2012_withAnnot)

sub_outliers_2002and2012_withAnnot <- subset(outliers_2002and2012_withAnnot, V5 > (WinStart - 50000) & V4 < (WinStart + 50000))
for_FileS4_2002and2012 <- sub_outliers_2002and2012_withAnnot[,c(1,2,7,15,16,18,20)]
colnames(for_FileS4_2002and2012) <- c("Scaf_ID", "Scaf", "Chr", "start", "stop", "strand", "gene")

write.table(for_File4_2002and2012 %>% distinct(), file = "Genes_50kb_outliers_2002&2012.txt") 


#2012-2017
outliers_2012and2017 <- rbind(hiCry1FST2, hiCry2FST2)
outliers_2012and2017_withIDs <- merge(outliers_2012and2017, scaf_names, by.x = "Scaf", by.y = "CHROM")

uniqScafs_2012and2017 <- as.character(unique(outliers_2012and2017_withIDs$Scaf_ID))
uniqScafs_2012and2017

all_genes_outliers_2012and2017 <- data.frame()

for (item in 1:length(uniqScafs_2012and2017)){
  sample <- uniqScafs_2012and2017[item]
  print(sample)
  Genes_2012and2017 <- ordered_gff3_genes[grepl(sample, ordered_gff3_genes$V1), ]
  all_genes_outliers_2012and2017 <- rbind(Genes_2012and2017,all_genes_outliers_2012and2017)
}

outliers_2012and2017_withAnnot <- merge(outliers_2012and2017_withIDs, all_genes_outliers_2012and2017, by.x = "Scaf_ID", by.y = "V1")
str(outliers_2012and2017_withAnnot)

sub_outliers_2012and2017_withAnnot <- subset(outliers_2012and2017_withAnnot, V5 > (WinStart - 50000) & V4 < (WinStart + 50000))
for_FileS4_2012and2017 <- sub_outliers_2012and2017_withAnnot[,c(1,2,7,15,16,18,20)]
colnames(for_FileS4_2012and2017) <- c("Scaf_ID", "Scaf", "Chr", "start", "stop", "strand", "gene")

write.table(for_FileS4_2012and2017 %>% distinct(), file = "Genes_50kb_outliers_2012and2017.txt") 

