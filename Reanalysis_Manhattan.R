#Here is the script I used to order and orient scaffolds and generate new Manhattan plots
#03162021 MF

#generating CMplots and getting scaffolds with outliers.

library(LEA);library(adegenet);library(vcfR);library(ape);library(CMplot); library(plyr)

setwd("~/Desktop/Hz_fieldColl_pop_gen/Reanalysis_PNAS")

#loading chromosomes
LG_ord_ZeaScafs <- read.csv("DEO9A1_scaffold_order.csv")
head(LG_ord_ZeaScafs)

Arm_ZeaScafs <- read.table("hzea.scaf.arm.agp")

#removing gaps
ZeaScafs_clean <- subset(Arm_ZeaScafs, V5 == "W")
head(ZeaScafs_clean)

#merging agp and LG files to get ordered arm scafs
Chroms <- data.frame(merge(ZeaScafs_clean, LG_ord_ZeaScafs, by.x = "V6", by.y = "scaf"))
head(Chroms)
str(Chroms)

ord_Chroms <- Chroms[order(Chroms[,11], Chroms[, 12], Chroms[, 2], Chroms[, 3]),]

#exporting for visual inspection and some minor adjustments for arm scaf order.
#Where an arm scaf was split toward the end, adjusted order to keep arm scafs together

#write.csv(ord_Chroms, file = "ord_Chroms_forEdit.csv", row.names = F)

#reading edited file back in
ord_Chroms_ed <- read.csv("ord_Chroms_forEdit.csv", header = T)
head(ord_Chroms_ed)

#sanity check that the loop works.  Pulled the first row and copied it to bottom row and made it part of chr31
ord_Chroms_ed_test <- read.csv("ord_Chroms_forEdit_test.csv", header = T)
head(ord_Chroms_ed_test)
tail(ord_Chroms_ed_test)#note first and last line are the same.

#final data after running loop and correcting
ord_Chroms_final <- read.csv("ord_Chroms_final.csv", header = T)

#Loop for checking for duplicate arm scaffolds
data <- data.frame(ord_Chroms_final)#change dataset when not testing. 
dups <- vector()

for (i in seq(c(1:31))){
  print(i)
  sub_chrom <- subset(data, chr == i)
  sub_els <- subset(data, chr != i)
  uniq_chrom <- unique(sub_chrom$V1)
  uniq_els <- unique(sub_els$V1)
  dups <- c(dups, intersect(uniq_chrom, uniq_els))
  
}

#getting duplicated scafs - for "final" data the counts table should be empty
uniq_dups <- data.frame(unique(dups))#this gives the list of duplicated armigera scafs in the linkage map.
names(uniq_dups)[1] <- "V1"

dups_fullData <- merge(uniq_dups, data, by = "V1")

counts <- ddply(dups_fullData, .(dups_fullData$V1, dups_fullData$chr), nrow)
names(counts) <- c("V1", "chr", "Freq")
print(counts)#Getting a table of the number of occurrences per chromosome of duplicates

#in most cases a clear choice - keeping the scaffold on the chr where it has the most hits to markers (Freq)
#in cases where there is not a clear choice - calling scafs a separate chromosome


#Use final clean dataset to order the scaffolds
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


#first 2002-2017 comparison
#40kb
wcFST_2002_2017_unfilt <- read.table("2002and2017_40kb_wcFST_all.smoothed", header = F)
names(wcFST_2002_2017_unfilt) <- c("Scaf", "WinStart", "WinStop", "NumSnps", "wcFST")
wcFST_2002_2017_unfilt$Scaf <- as.character(wcFST_2002_2017_unfilt$Scaf)

wcFST_2002_2017 <- subset(wcFST_2002_2017_unfilt, "NumSnps" > 10)

#merging the linkage map ordered scaffolds 

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


#2002-2012 comparison
#40kb
wcFST_2002_2012_unfilt <- read.table("2002and2012_40kb_wcFST_all.smoothed", header = F)
names(wcFST_2002_2012_unfilt) <- c("Scaf", "WinStart", "WinStop", "NumSnps", "wcFST")
wcFST_2002_2012_unfilt$Scaf <- as.character(wcFST_2002_2012_unfilt$Scaf)

wcFST_2002_2012 <- subset(wcFST_2002_2012_unfilt, "NumSnps" > 10)

#getting significance threshold based on ztransformation
wcFST_2002_2012$ztrans <- scale(wcFST_2002_2012$wcFST, center = TRUE, scale = TRUE) #ztransformation 
wcFST_2002_2012$SnpName <- seq(1:NROW(wcFST_2002_2012))

hist(wcFST_2002_2012$ztrans, breaks = 100)

hi_wcFST0 <- subset(wcFST_2002_2012, ztrans > 6) #subsetting out any window with an FST greater than 6 SDs from mean.
print(min(hi_wcFST0$wcFST))
print(mean(wcFST_2002_2012$wcFST))
print(median(wcFST_2002_2012$wcFST))

#merging ordered scafs with SNPs for plot
forPlot_2002_2012 <- data.frame(cbind(wcFST_2002_2012$SnpName,wcFST_2002_2012$Scaf, wcFST_2002_2012$WinStart, wcFST_2002_2012$wcFST))
names(forPlot_2002_2012) <- c("SnpName", "Scaf", "WinStart", "wcFST")

merged_2002_2012 <- merge(Full_Ord_Scafs, forPlot_2002_2012, by = "Scaf")
merged_2002_2012 <- merged_2002_2012[order(merged_2002_2012[,5], merged_2002_2012[,6], merged_2002_2012[,3], 
                                           merged_2002_2012[,1], merged_2002_2012[,8]),]

#adding artificial continuous window
with_artWin <- data.frame()

for (i in seq(c(1:max(merged_2002_2012$Chr)))){
  print(i)
  sub_chrom <- subset(merged_2002_2012, Chr == i)
  artWin <- seq(0, (nrow(sub_chrom)-1)*10000, by = 10000)
  comb <- data.frame(cbind(sub_chrom, artWin))
  with_artWin <- rbind(with_artWin, comb)
}

merged_2002_2012_forPlot <- with_artWin[,c(7,5,10,9)]
names(merged_2002_2012_forPlot) <- c("SnpName", "Chr", "WinStart", "wcFST")

#2002-2012 Manhattan plot function
CMplot(merged_2002_2012_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST0$wcFST)),
       chr.den.col="pink", file="jpg", memo="2002and2012_40kb_wcFST", dpi=300, threshold = min(hi_wcFST0$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hi_wcFST0$SnpName)


#2012-2017 comparison
#40kb
wcFST_2012_2017_unfilt <- read.table("2012and2017_40kb_wcFST_all.smoothed", header = F)
names(wcFST_2012_2017_unfilt) <- c("Scaf", "WinStart", "WinStop", "NumSnps", "wcFST")
wcFST_2012_2017_unfilt$Scaf <- as.character(wcFST_2012_2017_unfilt$Scaf)

wcFST_2012_2017 <- subset(wcFST_2012_2017_unfilt, "NumSnps" > 10)

#getting significance threshold based on ztransformation
wcFST_2012_2017$ztrans <- scale(wcFST_2012_2017$wcFST, center = TRUE, scale = TRUE) #ztransformation 
wcFST_2012_2017$SnpName <- seq(1:NROW(wcFST_2012_2017))

hist(wcFST_2012_2017$ztrans, breaks = 100)

hi_wcFST0 <- subset(wcFST_2012_2017, ztrans > 6) #subsetting out any window with an FST greater than 6 SDs from mean.
print(min(hi_wcFST0$wcFST))
print(mean(wcFST_2012_2017$wcFST))
print(median(wcFST_2012_2017$wcFST))

#merging ordered scafs with SNPs for plot
forPlot_2012_2017 <- data.frame(cbind(wcFST_2012_2017$SnpName,wcFST_2012_2017$Scaf, wcFST_2012_2017$WinStart, wcFST_2012_2017$wcFST))
names(forPlot_2012_2017) <- c("SnpName", "Scaf", "WinStart", "wcFST")

merged_2012_2017 <- merge(Full_Ord_Scafs, forPlot_2012_2017, by = "Scaf")
merged_2012_2017 <- merged_2012_2017[order(merged_2012_2017[,5], merged_2012_2017[,6], merged_2012_2017[,3], 
                                           merged_2012_2017[,1], merged_2012_2017[,8]),]

#adding artificial continuous window
with_artWin <- data.frame()

for (i in seq(c(1:max(merged_2012_2017$Chr)))){
  print(i)
  sub_chrom <- subset(merged_2012_2017, Chr == i)
  artWin <- seq(0, (nrow(sub_chrom)-1)*10000, by = 10000)
  comb <- data.frame(cbind(sub_chrom, artWin))
  with_artWin <- rbind(with_artWin, comb)
}

merged_2012_2017_forPlot <- with_artWin[,c(7,5,10,9)]
names(merged_2012_2017_forPlot) <- c("SnpName", "Chr", "WinStart", "wcFST")

#2012-2017 Manhattan plot function
CMplot(merged_2012_2017_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST0$wcFST)),
       chr.den.col="pink", file="jpg", memo="2012and2017_40kb_wcFST", dpi=300, threshold = min(hi_wcFST0$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hi_wcFST0$SnpName)
