#Script to generate new Manhattan plots with 5-40kb windows - All 2002-2017
#03222021 MF

#library(LEA);library(adegenet);library(vcfR);library(ape);
library(CMplot); library(plyr)

ord_Chroms_final <- read.csv("~/Desktop/Hz_fieldColl_pop_gen/Reanalysis_PNAS/ord_Chroms_final.csv", header = T)

scafs_by_ChrMark <- ord_Chroms_final[,c(2,11)]
names(scafs_by_ChrMark) <- c("Scaffold", "Chr")

#loop to add order
data_with_ordChr <- data.frame()

for (i in seq(c(1:max(scafs_by_ChrMark$Chr)))){
  print(i)
  sub_chrom <- subset(scafs_by_ChrMark, Chr == i)
  uniq_scafs <- unique(sub_chrom[c("Scaffold", "Chr")])
  uniq_scafs $scaf_ord <- cbind(seq(c(1:nrow(uniq_scafs ))))
  data_with_ordChr <- rbind(data_with_ordChr, uniq_scafs)
}

Full_Scafs <- merge(ZeaScafs_clean, data_with_ordChr, by.x = "V1", by.y = "Scaffold")
Full_Ord_Scafs <- Full_Scafs[order(Full_Scafs[,10], Full_Scafs[, 11], Full_Scafs[, 2]),]
Full_Ord_Scafs <- Full_Ord_Scafs[,c(-4,-5,-7,-8,-9)]
names(Full_Ord_Scafs) <- c("ArmScaf", "Start", "Stop", "Scaf", "Chr", "ArmScaf_Ord" )


##### 2002-2017 comparison - 40kb #####
wcFST_2002_2017_unfilt <- read.table("2002and2017_40kb_wcFST_all.smoothed", header = F)
names(wcFST_2002_2017_unfilt) <- c("Scaf", "WinStart", "WinStop", "NumSnps", "wcFST")
wcFST_2002_2017_unfilt$Scaf <- as.character(wcFST_2002_2017_unfilt$Scaf)

wcFST_2002_2017 <- subset(wcFST_2002_2017_unfilt, "NumSnps" > 10)

#getting significance threshold based on ztransformation
wcFST_2002_2017$ztrans <- scale(wcFST_2002_2017$wcFST, center = TRUE, scale = TRUE) #ztransformation 
wcFST_2002_2017$SnpName <- seq(1:NROW(wcFST_2002_2017))

hist(wcFST_2002_2017$ztrans, breaks = 100)

hi_wcFST0 <- subset(wcFST_2002_2017, ztrans > 6) #subsetting out any window with an FST greater than 6 SDs from mean.
print(min(hi_wcFST0$wcFST))
print(mean(wcFST_2002_2017$wcFST))
print(median(wcFST_2002_2017$wcFST))

#merging ordered scafs with SNPs for plot
forPlot_2002_2017 <- data.frame(cbind(wcFST_2002_2017$SnpName,wcFST_2002_2017$Scaf, wcFST_2002_2017$WinStart, wcFST_2002_2017$wcFST))
names(forPlot_2002_2017) <- c("SnpName", "Scaf", "WinStart", "wcFST")

merged_2002_2017 <- merge(Full_Ord_Scafs, forPlot_2002_2017, by = "Scaf")
merged_2002_2017 <- merged_2002_2017[order(merged_2002_2017 [,5], merged_2002_2017 [,6], merged_2002_2017 [,3], 
                                           merged_2002_2017 [,1], merged_2002_2017 [,8]),]

#adding artificial continuous window
with_artWin <- data.frame()

for (i in seq(c(1:max(merged_2002_2017$Chr)))){
  print(i)
  sub_chrom <- subset(merged_2002_2017, Chr == i)
  artWin <- seq(0, (nrow(sub_chrom)-1)*10000, by = 10000)
  comb <- data.frame(cbind(sub_chrom, artWin))
  with_artWin <- rbind(with_artWin, comb)
}

merged_2002_2017_forPlot <- with_artWin[,c(7,5,10,9)]
names(merged_2002_2017_forPlot) <- c("SnpName", "Chr", "WinStart", "wcFST")

#2002-2017 Manhattan plot function
CMplot(merged_2002_2017_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST0$wcFST)),
       chr.den.col="pink", file="jpg", memo="2002and2017_40kb_wcFST", dpi=300, threshold = min(hi_wcFST0$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hi_wcFST0$SnpName)



##### 2002-2017 comparison - 20kb #####
wcFST_2002_2017_unfilt <- read.table("2002and2017_20kb_wcFST_all.smoothed", header = F)
names(wcFST_2002_2017_unfilt) <- c("Scaf", "WinStart", "WinStop", "NumSnps", "wcFST")
wcFST_2002_2017_unfilt$Scaf <- as.character(wcFST_2002_2017_unfilt$Scaf)

wcFST_2002_2017 <- subset(wcFST_2002_2017_unfilt, "NumSnps" > 10)

#getting significance threshold based on ztransformation
wcFST_2002_2017$ztrans <- scale(wcFST_2002_2017$wcFST, center = TRUE, scale = TRUE) #ztransformation 
wcFST_2002_2017$SnpName <- seq(1:NROW(wcFST_2002_2017))

hist(wcFST_2002_2017$ztrans, breaks = 100)

hi_wcFST0 <- subset(wcFST_2002_2017, ztrans > 6) #subsetting out any window with an FST greater than 6 SDs from mean.
print(min(hi_wcFST0$wcFST))
print(mean(wcFST_2002_2017$wcFST))
print(median(wcFST_2002_2017$wcFST))

#merging ordered scafs with SNPs for plot
forPlot_2002_2017 <- data.frame(cbind(wcFST_2002_2017$SnpName,wcFST_2002_2017$Scaf, wcFST_2002_2017$WinStart, wcFST_2002_2017$wcFST))
names(forPlot_2002_2017) <- c("SnpName", "Scaf", "WinStart", "wcFST")

merged_2002_2017 <- merge(Full_Ord_Scafs, forPlot_2002_2017, by = "Scaf")
merged_2002_2017 <- merged_2002_2017[order(merged_2002_2017 [,5], merged_2002_2017 [,6], merged_2002_2017 [,3], 
                                           merged_2002_2017 [,1], merged_2002_2017 [,8]),]

#adding artificial continuous window
with_artWin <- data.frame()

for (i in seq(c(1:max(merged_2002_2017$Chr)))){
  print(i)
  sub_chrom <- subset(merged_2002_2017, Chr == i)
  artWin <- seq(0, (nrow(sub_chrom)-1)*10000, by = 10000)
  comb <- data.frame(cbind(sub_chrom, artWin))
  with_artWin <- rbind(with_artWin, comb)
}

merged_2002_2017_forPlot <- with_artWin[,c(7,5,10,9)]
names(merged_2002_2017_forPlot) <- c("SnpName", "Chr", "WinStart", "wcFST")

#2002-2017 Manhattan plot function
CMplot(merged_2002_2017_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST0$wcFST)),
       chr.den.col="pink", file="jpg", memo="2002and2017_20kb_wcFST", dpi=300, threshold = min(hi_wcFST0$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hi_wcFST0$SnpName)


##### 2002-2017 comparison - 10kb #####
wcFST_2002_2017_unfilt <- read.table("2002and2017_10kb_wcFST_all.smoothed", header = F)
names(wcFST_2002_2017_unfilt) <- c("Scaf", "WinStart", "WinStop", "NumSnps", "wcFST")
wcFST_2002_2017_unfilt$Scaf <- as.character(wcFST_2002_2017_unfilt$Scaf)

wcFST_2002_2017 <- subset(wcFST_2002_2017_unfilt, "NumSnps" > 10)

#getting significance threshold based on ztransformation
wcFST_2002_2017$ztrans <- scale(wcFST_2002_2017$wcFST, center = TRUE, scale = TRUE) #ztransformation 
wcFST_2002_2017$SnpName <- seq(1:NROW(wcFST_2002_2017))

hist(wcFST_2002_2017$ztrans, breaks = 100)

hi_wcFST0 <- subset(wcFST_2002_2017, ztrans > 6) #subsetting out any window with an FST greater than 6 SDs from mean.
print(min(hi_wcFST0$wcFST))
print(mean(wcFST_2002_2017$wcFST))
print(median(wcFST_2002_2017$wcFST))

#merging ordered scafs with SNPs for plot
forPlot_2002_2017 <- data.frame(cbind(wcFST_2002_2017$SnpName,wcFST_2002_2017$Scaf, wcFST_2002_2017$WinStart, wcFST_2002_2017$wcFST))
names(forPlot_2002_2017) <- c("SnpName", "Scaf", "WinStart", "wcFST")

merged_2002_2017 <- merge(Full_Ord_Scafs, forPlot_2002_2017, by = "Scaf")
merged_2002_2017 <- merged_2002_2017[order(merged_2002_2017 [,5], merged_2002_2017 [,6], merged_2002_2017 [,3], 
                                           merged_2002_2017 [,1], merged_2002_2017 [,8]),]

#adding artificial continuous window
with_artWin <- data.frame()

for (i in seq(c(1:max(merged_2002_2017$Chr)))){
  print(i)
  sub_chrom <- subset(merged_2002_2017, Chr == i)
  artWin <- seq(0, (nrow(sub_chrom)-1)*10000, by = 10000)
  comb <- data.frame(cbind(sub_chrom, artWin))
  with_artWin <- rbind(with_artWin, comb)
}

merged_2002_2017_forPlot <- with_artWin[,c(7,5,10,9)]
names(merged_2002_2017_forPlot) <- c("SnpName", "Chr", "WinStart", "wcFST")

#2002-2017 Manhattan plot function
CMplot(merged_2002_2017_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST0$wcFST)),
       chr.den.col="pink", file="jpg", memo="2002and2017_10kb_wcFST", dpi=300, threshold = min(hi_wcFST0$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hi_wcFST0$SnpName)


##### 2002-2017 comparison - 5kb #####
wcFST_2002_2017_unfilt <- read.table("2002and2017_5kb_wcFST_all.smoothed", header = F)
names(wcFST_2002_2017_unfilt) <- c("Scaf", "WinStart", "WinStop", "NumSnps", "wcFST")
wcFST_2002_2017_unfilt$Scaf <- as.character(wcFST_2002_2017_unfilt$Scaf)

wcFST_2002_2017 <- subset(wcFST_2002_2017_unfilt, "NumSnps" > 10)

#getting significance threshold based on ztransformation
wcFST_2002_2017$ztrans <- scale(wcFST_2002_2017$wcFST, center = TRUE, scale = TRUE) #ztransformation 
wcFST_2002_2017$SnpName <- seq(1:NROW(wcFST_2002_2017))

hist(wcFST_2002_2017$ztrans, breaks = 100)

hi_wcFST0 <- subset(wcFST_2002_2017, ztrans > 6) #subsetting out any window with an FST greater than 6 SDs from mean.
print(min(hi_wcFST0$wcFST))
print(mean(wcFST_2002_2017$wcFST))
print(median(wcFST_2002_2017$wcFST))

#merging ordered scafs with SNPs for plot
forPlot_2002_2017 <- data.frame(cbind(wcFST_2002_2017$SnpName,wcFST_2002_2017$Scaf, wcFST_2002_2017$WinStart, wcFST_2002_2017$wcFST))
names(forPlot_2002_2017) <- c("SnpName", "Scaf", "WinStart", "wcFST")

merged_2002_2017 <- merge(Full_Ord_Scafs, forPlot_2002_2017, by = "Scaf")
merged_2002_2017 <- merged_2002_2017[order(merged_2002_2017 [,5], merged_2002_2017 [,6], merged_2002_2017 [,3], 
                                           merged_2002_2017 [,1], merged_2002_2017 [,8]),]

#adding artificial continuous window
with_artWin <- data.frame()

for (i in seq(c(1:max(merged_2002_2017$Chr)))){
  print(i)
  sub_chrom <- subset(merged_2002_2017, Chr == i)
  artWin <- seq(0, (nrow(sub_chrom)-1)*10000, by = 10000)
  comb <- data.frame(cbind(sub_chrom, artWin))
  with_artWin <- rbind(with_artWin, comb)
}

merged_2002_2017_forPlot <- with_artWin[,c(7,5,10,9)]
names(merged_2002_2017_forPlot) <- c("SnpName", "Chr", "WinStart", "wcFST")

#2002-2017 Manhattan plot function
CMplot(merged_2002_2017_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST0$wcFST)),
       chr.den.col="pink", file="jpg", memo="2002and2017_5kb_wcFST", dpi=300, threshold = min(hi_wcFST0$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hi_wcFST0$SnpName)