#Script to generate new Manhattan plots with 5-40kb windows - All 2002-2017
#03222021 MF

#library(LEA);library(adegenet);library(vcfR);library(ape);
library(CMplot);library(dplyr)

setwd("~/Desktop/Hz_fieldColl_pop_gen/Reanalysis_PNAS/")

#loading zea superscaffolds
Full_Ord_Scafs <- read.table(file = "Hzea_superScaf_genome.txt", header = T)

##### 2002-2017 comparison - 40kb with 10kb step #####
wcFST_2002_2017_unfilt <- read.table("2002and2017_40kb_wcFST_all.smoothed", header = F)
names(wcFST_2002_2017_unfilt) <- c("Scaf", "WinStart", "WinStop", "NumSnps", "wcFST")
wcFST_2002_2017_unfilt$Scaf <- as.character(wcFST_2002_2017_unfilt$Scaf)

wcFST_2002_2017 <- subset(wcFST_2002_2017_unfilt, NumSnps > 10)

#getting significance threshold based on ztransformation
ztrans <- as.vector(scale(wcFST_2002_2017$wcFST, center = TRUE, scale = TRUE)) #ztransformation
wcFST_2002_2017 <- data.frame(cbind(wcFST_2002_2017, ztrans))
wcFST_2002_2017$SnpName <- seq(1:NROW(wcFST_2002_2017))

hist(wcFST_2002_2017$ztrans, breaks = 100)

hi_wcFST0 <- subset(wcFST_2002_2017, ztrans > 6) #subsetting out any window with an FST greater than 6 SDs from mean.
print(min(hi_wcFST0$wcFST))
print(mean(wcFST_2002_2017$wcFST))
print(median(wcFST_2002_2017$wcFST))

#merging ordered scafs with SNPs for plot
forPlot_2002_2017 <- data.frame(cbind(wcFST_2002_2017$SnpName,wcFST_2002_2017$Scaf, wcFST_2002_2017$WinStart, wcFST_2002_2017$wcFST))
names(forPlot_2002_2017) <- c("SnpName", "Scaf", "WinStart", "wcFST")
forPlot_2002_2017$wcFST <- as.numeric(as.character(forPlot_2002_2017$wcFST))
forPlot_2002_2017$WinStart <- as.numeric(as.character(forPlot_2002_2017$WinStart))

merged_2002_2017 <- merge(Full_Ord_Scafs, forPlot_2002_2017, by = "Scaf")
merged_2002_2017 <- merged_2002_2017[order(merged_2002_2017 [,5], merged_2002_2017 [,6], merged_2002_2017 [,3], 
                                           merged_2002_2017 [,1], merged_2002_2017 [,8]),]

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

merged_2002_2017_forPlot <- with_artWin[,c(7,5,11,9)]
names(merged_2002_2017_forPlot) <- c("SnpName", "Chr", "WinStart", "wcFST")


#2002-2017 Manhattan plot function
CMplot(merged_2002_2017_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST0$wcFST)),
       chr.den.col="pink", file="jpg", memo="2002and2017_40kb_wcFST", dpi=300, threshold = min(hi_wcFST0$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hi_wcFST0$SnpName)



##### 2002-2017 comparison - 20kb with 5kb step #####
wcFST_2002_2017_unfilt <- read.table("2002and2017_20kb_wcFST_all.smoothed", header = F)
names(wcFST_2002_2017_unfilt) <- c("Scaf", "WinStart", "WinStop", "NumSnps", "wcFST")
wcFST_2002_2017_unfilt$Scaf <- as.character(wcFST_2002_2017_unfilt$Scaf)

wcFST_2002_2017 <- subset(wcFST_2002_2017_unfilt, NumSnps > 10)

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
forPlot_2002_2017$wcFST <- as.numeric(as.character(forPlot_2002_2017$wcFST))
forPlot_2002_2017$WinStart <- as.numeric(as.character(forPlot_2002_2017$WinStart))

merged_2002_2017 <- merge(Full_Ord_Scafs, forPlot_2002_2017, by = "Scaf")
merged_2002_2017 <- merged_2002_2017[order(merged_2002_2017 [,5], merged_2002_2017 [,6], merged_2002_2017 [,3], 
                                           merged_2002_2017 [,1], merged_2002_2017 [,8]),]

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

merged_2002_2017_forPlot <- with_artWin[,c(7,5,11,9)]
names(merged_2002_2017_forPlot) <- c("SnpName", "Chr", "WinStart", "wcFST")

#2002-2017 Manhattan plot function
CMplot(merged_2002_2017_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST0$wcFST)),
       chr.den.col="pink", file="jpg", memo="2002and2017_20kb_wcFST", dpi=300, threshold = min(hi_wcFST0$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hi_wcFST0$SnpName)


##### 2002-2017 comparison - 10kb with 1kb step #####
wcFST_2002_2017_unfilt <- read.table("2002and2017_10kb_wcFST_all.smoothed", header = F)
names(wcFST_2002_2017_unfilt) <- c("Scaf", "WinStart", "WinStop", "NumSnps", "wcFST")
wcFST_2002_2017_unfilt$Scaf <- as.character(wcFST_2002_2017_unfilt$Scaf)

wcFST_2002_2017 <- subset(wcFST_2002_2017_unfilt, NumSnps > 10)

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
forPlot_2002_2017$wcFST <- as.numeric(as.character(forPlot_2002_2017$wcFST))
forPlot_2002_2017$WinStart <- as.numeric(as.character(forPlot_2002_2017$WinStart))

merged_2002_2017 <- merge(Full_Ord_Scafs, forPlot_2002_2017, by = "Scaf")
merged_2002_2017 <- merged_2002_2017[order(merged_2002_2017 [,5], merged_2002_2017 [,6], merged_2002_2017 [,3], 
                                           merged_2002_2017 [,1], merged_2002_2017 [,8]),]

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

merged_2002_2017_forPlot <- with_artWin[,c(7,5,11,9)]
names(merged_2002_2017_forPlot) <- c("SnpName", "Chr", "WinStart", "wcFST")

#2002-2017 Manhattan plot function
CMplot(merged_2002_2017_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST0$wcFST)),
       chr.den.col="pink", file="jpg", memo="2002and2017_10kb_wcFST", dpi=300, threshold = min(hi_wcFST0$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hi_wcFST0$SnpName)


##### 2002-2017 comparison - 5kb with 1kb step #####
wcFST_2002_2017_unfilt <- read.table("2002and2017_5kb_wcFST_all.smoothed", header = F)
names(wcFST_2002_2017_unfilt) <- c("Scaf", "WinStart", "WinStop", "NumSnps", "wcFST")
wcFST_2002_2017_unfilt$Scaf <- as.character(wcFST_2002_2017_unfilt$Scaf)

wcFST_2002_2017 <- subset(wcFST_2002_2017_unfilt, NumSnps > 10)

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
forPlot_2002_2017$wcFST <- as.numeric(as.character(forPlot_2002_2017$wcFST))
forPlot_2002_2017$WinStart <- as.numeric(as.character(forPlot_2002_2017$WinStart))

merged_2002_2017 <- merge(Full_Ord_Scafs, forPlot_2002_2017, by = "Scaf")
merged_2002_2017 <- merged_2002_2017[order(merged_2002_2017 [,5], merged_2002_2017 [,6], merged_2002_2017 [,3], 
                                           merged_2002_2017 [,1], merged_2002_2017 [,8]),]

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

merged_2002_2017_forPlot <- with_artWin[,c(7,5,11,9)]
names(merged_2002_2017_forPlot) <- c("SnpName", "Chr", "WinStart", "wcFST")

#2002-2017 Manhattan plot function
CMplot(merged_2002_2017_forPlot, plot.type="m", col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST0$wcFST)),
       chr.den.col="pink", file="jpg", memo="2002and2017_5kb_wcFST", dpi=300, threshold = min(hi_wcFST0$wcFST), LOG10 = F, ylab = "FST", xlab = "",
       highlight = hi_wcFST0$SnpName)