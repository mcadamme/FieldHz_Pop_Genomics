#Here is the script I used to generate new Manhattan plots with 10kb window sizes and add QTL
#03162021 MF

#generating CMplots and getting scaffolds with outliers.

library(CMplot)#; library(plyr)

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
print(min(hi_wcFST0$wcFST))
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
hi_wcFST0 <- subset(merged_2002_2012, ztrans > 6) #subsetting out any window with an FST greater than 6 SDs from mean.
print(min(hi_wcFST0$wcFST))
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
hiCry1FST <- subset(Cry1FST, ztrans > 6)

write.table(hiCry1FST, "hiCry1FST2002to2012.txt", row.names = F)

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
hiCry2FST <- subset(Cry2FST, ztrans > 6)

write.table(hiCry2FST, "hiCry2FST2002to2012.txt", row.names = F)

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
CMplot(merged_2002_2012_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST0$wcFST)),
       chr.den.col="pink", file="jpg", memo="2002and2012_10kb_wcFST", dpi=300, threshold = min(hi_wcFST0$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hi_wcFST0$SnpName)

#sig Cry1
CMplot(merged_2002_2012_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST0$wcFST)),
       chr.den.col="pink", file="jpg", memo="2002and2012_10kb_cry1", dpi=300, threshold = min(hi_wcFST0$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hiCry1FST$SnpName)

#sig Cry2
CMplot(merged_2002_2012_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST0$wcFST)),
       chr.den.col="pink", file="jpg", memo="2002and2012_10kb_cry2", dpi=300, threshold = min(hi_wcFST0$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hiCry2FST$SnpName)


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
hi_wcFST0 <- subset(merged_2012_2017, ztrans > 6) #subsetting out any window with an FST greater than 6 SDs from mean.
print(min(hi_wcFST0$wcFST))
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
hiCry1FST <- subset(Cry1FST, ztrans > 6)

write.table(hiCry1FST, "hiCry1FST2012to2017.txt", row.names = F)

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
hiCry2FST <- subset(Cry2FST, ztrans > 6)

write.table(hiCry2FST, "hiCry2FST2012to2017.txt", row.names = F)

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
CMplot(merged_2012_2017_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST0$wcFST)),
       chr.den.col="pink", file="jpg", memo="2012and2017_10kb_wcFST", dpi=300, threshold = min(hi_wcFST0$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hi_wcFST0$SnpName)

#sig Cry1
CMplot(merged_2012_2017_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST0$wcFST)),
       chr.den.col="pink", file="jpg", memo="2012and2017_10kb_cry1", dpi=300, threshold = min(hi_wcFST0$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hiCry1FST$SnpName)

#sig Cry2
CMplot(merged_2012_2017_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST0$wcFST)),
       chr.den.col="pink", file="jpg", memo="2012and2017_10kb_cry2", dpi=300, threshold = min(hi_wcFST0$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hiCry2FST$SnpName)

#What are the scaffolds not represented in this final plot?

scaf_names <- read.table("/media/megan/New Volume1/Hz_PopGen_ddRAD_demult/Bowtie_genome_alignments/Hzea_genome/scaffold_and_contig_names.txt",
                         header = T)

all_scafs <- unique(scaf_names$CHROM)
plotted_scafs <- unique(Full_Ord_Scafs$Scaf)
diff <- setdiff(all_scafs, plotted_scafs)
diff_df <- data.frame(diff)

#no candidate genes from Table 1 on these scaffolds - checked by hand.

intersect(diff_df$diff, unique(hi_wcFST0$Scaf))#as of 3/22 this should be 0.
#added "KZ116441.1" "KZ116446.1" "KZ117435.1" to ord_Chroms_final on 3/22/21, as chromsomes 43-45.  
#These could not be mapped to chromosomes, even by hand.  No known genes on them. Very short.

