#Script for super-scaffolding zea genome using .agp file from armigera alignment and the linkage map for DE109.
#03232021 MF

library(plyr)

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

#final data after running loop and correcting plus added scaffold containing candidate genes that were missing.

#in most cases a clear choice - kept scaffold on the chr where it has the most hits to markers (Freq)
#in cases where there is not a clear choice - calling scafs a separate chromosome

ord_Chroms_final <- read.csv("ord_Chroms_final.csv", header = T)
str(ord_Chroms_final)
table(ord_Chroms_final$chr)#shows 46 "chromosomes" as expected


#Loop for checking for duplicate arm scaffolds
data <- data.frame(ord_Chroms_final)#change dataset when not testing. 
dups <- vector()

for (i in seq(c(1:46))){
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


#Use final clean dataset to order the scaffolds
scafs_by_ChrMark <- ord_Chroms_final[,c(2,11)]
names(scafs_by_ChrMark) <- c("Scaffold", "Chr")
table(scafs_by_ChrMark$Chr)#sanity check - shows 46 chr

#loop to add order
data_with_ordChr <- data.frame()

for (i in seq(c(1:max(scafs_by_ChrMark$Chr)))){
  print(i)
  sub_chrom <- subset(scafs_by_ChrMark, Chr == i)
  uniq_scafs <- unique(sub_chrom[c("Scaffold", "Chr")])
  uniq_scafs$scaf_ord <- cbind(seq(c(1:nrow(uniq_scafs ))))
  data_with_ordChr <- rbind(data_with_ordChr, uniq_scafs)
}

str(data_with_ordChr)
table(data_with_ordChr$Chr)

Full_Scafs <- merge(ZeaScafs_clean, data_with_ordChr, by.x = "V1", by.y = "Scaffold")
Full_Ord_Scafs <- Full_Scafs[order(Full_Scafs[,10], Full_Scafs[, 11], Full_Scafs[, 2]),]
Full_Ord_Scafs <- Full_Ord_Scafs[,c(-4,-5,-7,-8,-9)]
names(Full_Ord_Scafs) <- c("ArmScaf", "Start", "Stop", "Scaf", "Chr", "ArmScaf_Ord" )
table(Full_Ord_Scafs$Chr)#sanity check - shows 46 chr
write.table(Full_Ord_Scafs, file = "Hzea_superScaf_genome.txt", row.names = F)
