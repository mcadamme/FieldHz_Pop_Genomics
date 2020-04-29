#generating CMplots and getting scaffolds with outliers.

library(LEA);library(adegenet);library(vcfR);library(CMplot)


setwd("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1")


#first 2002-2017 comparison
#40kb
wcFST_2002_2017_unfilt <- read.table("2002and2017_40kb_wcFST_all.smoothed", header = F)

wcFST_2002_2017 <- subset(wcFST_2002_2017_unfilt, V4 > 10)

#getting significance threshold based on ztransformation
wcFST_2002_2017$ztrans <- scale(wcFST_2002_2017$V5, center = TRUE, scale = TRUE) #ztransformation 
wcFST_2002_2017$SnpName <- seq(1:NROW(wcFST_2002_2017))

hist(wcFST_2002_2017$ztrans, breaks = 100)

hi_wcFST0 <- subset(wcFST_2002_2017, ztrans > 6) #subsetting out any window with an FST greater than 6 SDs from mean.
print(min(hi_wcFST0$V5))
print(mean(wcFST_2002_2017$V5))
print(median(wcFST_2002_2017$V5))

forPlot_2002_2017 <- data.frame(cbind(wcFST_2002_2017$SnpName,wcFST_2002_2017$V1, wcFST_2002_2017$V2, wcFST_2002_2017$V5))

#2002-2017 Circular Manhattan plot function
##IMPORTANT - each CMplot has to be run and fig must be renamed one at a time.

CMplot(forPlot_2002_2017, plot.type="c", r=1.6, cir.legend=TRUE, col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST0$V5)),
       outward=TRUE, cir.legend.col="black", cir.chr.h=.2 ,chr.den.col="pink", file="jpg",
       memo="2002and2017_40kb_wcFST", dpi=300, cex.lab = 4, threshold = min(hi_wcFST0$V5), LOG10 = F, chr.labels=seq(1,2958), 
       highlight = hi_wcFST0$SnpName)


#20kb
wcFST_2002_2017 <- read.table("2002and2017_20kb_wcFST_all.smoothed", header = F)

#getting significance threshold based on ztransformation
wcFST_2002_2017$ztrans <- scale(wcFST_2002_2017$V5, center = TRUE, scale = TRUE) #ztransformation 
wcFST_2002_2017$SnpName <- seq(1:NROW(wcFST_2002_2017))

hist(wcFST_2002_2017$ztrans, breaks = 100)

hi_wcFST1 <- subset(wcFST_2002_2017, ztrans > 6) #subsetting out any window with an FST greater than 6 SDs from mean.
print(min(hi_wcFST1$V5))


forPlot_2002_2017 <- data.frame(cbind(wcFST_2002_2017$SnpName,wcFST_2002_2017$V1, wcFST_2002_2017$V2, wcFST_2002_2017$V5))

#2002-2017 Circular Manhattan plot function
##IMPORTANT - each CMplot has to be run and fig must be renamed one at a time.

CMplot(forPlot_2002_2017, plot.type="c", r=1.6, cir.legend=TRUE, col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST1$V5)),
       outward=TRUE, cir.legend.col="black", cir.chr.h=.2 ,chr.den.col="pink", file="jpg",
       memo="2002and2017_20kb_wcFST", dpi=300, cex.lab = 4, threshold = min(hi_wcFST1$V5), LOG10 = F, chr.labels=seq(1,2958), 
       highlight = hi_wcFST1$SnpName)

#10kb
wcFST_2002_2017 <- read.table("2002and2017_10kb_wcFST_all.smoothed", header = F)

#getting significance threshold based on ztransformation
wcFST_2002_2017$ztrans <- scale(wcFST_2002_2017$V5, center = TRUE, scale = TRUE) #ztransformation 
wcFST_2002_2017$SnpName <- seq(1:NROW(wcFST_2002_2017))

hist(wcFST_2002_2017$ztrans, breaks = 100)

hi_wcFST2 <- subset(wcFST_2002_2017, ztrans > 6) #subsetting out any window with an FST greater than 6 SDs from mean.
print(min(hi_wcFST2$V5))


forPlot_2002_2017 <- data.frame(cbind(wcFST_2002_2017$SnpName,wcFST_2002_2017$V1, wcFST_2002_2017$V2, wcFST_2002_2017$V5))

#2002-2017 Circular Manhattan plot function
CMplot(forPlot_2002_2017, plot.type="c", r=1.6, cir.legend=TRUE, col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST2$V5)),
       outward=TRUE, cir.legend.col="black", cir.chr.h=.2 ,chr.den.col="pink", file="jpg",
       memo="2002and2017_10kb_wcFST", dpi=300, cex.lab = 4, threshold = min(hi_wcFST2$V5), LOG10 = F, chr.labels=seq(1,2958), 
       highlight = hi_wcFST2$SnpName)


#5kb
wcFST_2002_2017 <- read.table("2002and2017_5kb_wcFST_all.smoothed", header = F)

#getting significance threshold based on ztransformation
wcFST_2002_2017$ztrans <- scale(wcFST_2002_2017$V5, center = TRUE, scale = TRUE) #ztransformation 
wcFST_2002_2017$SnpName <- seq(1:NROW(wcFST_2002_2017))

hist(wcFST_2002_2017$ztrans, breaks = 100)

hi_wcFST3 <- subset(wcFST_2002_2017, ztrans > 6) #subsetting out any window with an FST greater than 6 SDs from mean.
print(min(hi_wcFST3$V5))


forPlot_2002_2017 <- data.frame(cbind(wcFST_2002_2017$SnpName,wcFST_2002_2017$V1, wcFST_2002_2017$V2, wcFST_2002_2017$V5))

#2002-2017 Circular Manhattan plot function
CMplot(forPlot_2002_2017, plot.type="c", r=1.6, cir.legend=TRUE, col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST3$V5)),
       outward=TRUE, cir.legend.col="black", cir.chr.h=.2 ,chr.den.col="pink", file="jpg",
       memo="2002and2017_5kb_wcFST", dpi=300, cex.lab = 4, threshold = min(hi_wcFST3$V5), LOG10 = F, chr.labels=seq(1,2958), 
       highlight = hi_wcFST3$SnpName)



#2002-2012 comparison
#40kb
wcFST_2002_2012_unfilt <- read.table("2002and2012_40kb_wcFST_all.smoothed", header = F)
wcFST_2002_2012 <- subset(wcFST_2002_2012_unfilt, V4 > 10)

#getting significance threshold based on ztransformation
wcFST_2002_2012$ztrans <- scale(wcFST_2002_2012$V5, center = TRUE, scale = TRUE) #ztransformation 
wcFST_2002_2012$SnpName <- seq(1:NROW(wcFST_2002_2012))

hist(wcFST_2002_2012$ztrans, breaks = 100)

hi_wcFST4 <- subset(wcFST_2002_2012, ztrans > 6) #subsetting out any window with an FST greater than 6 SDs from mean.
print(min(hi_wcFST4$V5))


forPlot_2002_2012 <- data.frame(cbind(wcFST_2002_2012$SnpName,wcFST_2002_2012$V1, wcFST_2002_2012$V2, wcFST_2002_2012$V5))

#2002-2012 Circular Manhattan plot function
CMplot(forPlot_2002_2012, plot.type="c", r=1.6, cir.legend=TRUE, col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,0.2),
       outward=TRUE, cir.legend.col="black", cir.chr.h=.2 ,chr.den.col="pink", file="jpg",
       memo="2002and2012_40kb_wcFST", dpi=300, cex.lab = 4, threshold = min(hi_wcFST4$V5), LOG10 = F, chr.labels=seq(1,2958), highlight = hi_wcFST4$SnpName)


#2012-2017 comparison
#40kb
wcFST_2012_2017_unfilt <- read.table("2012and2017_40kb_wcFST_all.smoothed", header = F)
wcFST_2012_2017 <- subset(wcFST_2012_2017_unfilt, V4 > 10)

#getting significance threshold based on ztransformation
wcFST_2012_2017$ztrans <- scale(wcFST_2012_2017$V5, center = TRUE, scale = TRUE) #ztransformation 
wcFST_2012_2017$SnpName <- seq(1:NROW(wcFST_2012_2017))

hist(wcFST_2012_2017$ztrans, breaks = 100)

hi_wcFST5 <- subset(wcFST_2012_2017, ztrans > 6) #subsetting out any window with an FST greater than 6 SDs from mean.
print(min(hi_wcFST5$V5))


forPlot_2012_2017 <- data.frame(cbind(wcFST_2012_2017$SnpName,wcFST_2012_2017$V1, wcFST_2012_2017$V2, wcFST_2012_2017$V5))

#2002-2012 Circular Manhattan plot function
CMplot(forPlot_2012_2017, plot.type="c", r=1.6, cir.legend=TRUE, col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,0.2),
       outward=TRUE, cir.legend.col="black", cir.chr.h=.2 ,chr.den.col="pink", file="jpg",
       memo="2012and2017_40kb_wcFST", dpi=300, cex.lab = 4, threshold = min(hi_wcFST5$V5), LOG10 = F, chr.labels=seq(1,2958), highlight = hi_wcFST5$SnpName)



#Loading scafname and gff files for getting genes nearby outliers

#loading identifier file
scaf_names <- read.table("/media/megan/New Volume1/Hz_PopGen_ddRAD_demult/Bowtie_genome_alignments/Hzea_genome/scaffold_and_contig_names.txt",
                         header = T)



#Now pulling out the parts of the gff3 file that are important
gff3 <- read.table("/media/megan/New Volume1/Hz_PopGen_ddRAD_demult/Bowtie_genome_alignments/Hzea_genome/HzOGS2-15205-fixed_note-added.gff3", sep="\t", stringsAsFactors = F)
head(gff3)


#comparisons for pairs of years - all 40kb windows
outliers_2002and2017_withIDs <- merge(hi_wcFST0, scaf_names, by.x = "V1", by.y = "CHROM")
write.table(outliers_2002and2017_withIDs, file = "2002and2017_outlierScaffolds.txt")
outliers_2002and2017_withAnnot <- merge(outliers_2002and2017_withIDs, gff3, by.x = "Scaf_ID", by.y = "V1")
sub_outliers_2002and2017_withAnnot <- subset(outliers_2002and2017_withAnnot, V3.y == "gene" & V5.y > (V2.x -50000) & V4.y < (V2.x + 50000))
write.table(sub_outliers_2002and2017_withAnnot, file = "Genes_50kb_outliers_2002&2017.txt")

outliers_2002and2012_withIDs <- merge(hi_wcFST4, scaf_names, by.x = "V1", by.y = "CHROM")
write.table(outliers_2002and2012_withIDs, file = "2002and2012_outlierScaffolds.txt")
outliers_2002and2012_withAnnot <- merge(outliers_2002and2012_withIDs, gff3, by.x = "Scaf_ID", by.y = "V1")
sub_outliers_2002and2012_withAnnot <- subset(outliers_2002and2012_withAnnot, V3.y == "gene" & V5.y > (V2.x -50000) & V4.y < (V2.x + 50000))
write.table(sub_outliers_2002and2012_withAnnot, file = "Genes_50kb_outliers_2002&2012.txt")

outliers_2012and2017_withIDs <- merge(hi_wcFST5, scaf_names, by.x = "V1", by.y = "CHROM")
write.table(outliers_2012and2017_withIDs, file = "2012and2017_outlierScaffolds.txt")
outliers_2012and2017_withAnnot <- merge(outliers_2012and2017_withIDs, gff3, by.x = "Scaf_ID", by.y = "V1")
sub_outliers_2012and2017_withAnnot <- subset(outliers_2012and2017_withAnnot, V3.y == "gene" & V5.y > (V2.x -50000) & V4.y < (V2.x + 50000))
write.table(sub_outliers_2012and2017_withAnnot, file = "Genes_50kb_outliers_2012&2017.txt")

#shared outliers
uniq_scafs_2002and2017 <- unique(sub_outliers_2002and2017_withAnnot$Scaf_ID)
uniq_scafs_2002and2017

uniq_scafs_2002and2012 <-unique(sub_outliers_2002and2012_withAnnot$Scaf_ID)
uniq_scafs_2002and2012

uniq_scafs_2012and2017 <-unique(sub_outliers_2012and2017_withAnnot$Scaf_ID)
uniq_scafs_2012and2017

#####Bt and Non-Bt Samples
setwd("/media/megan/New Volume1/Hzea_WGRS_Bowtie2_output/BtandNonBt_alignmentFiles/WGRS_mpileupANDvcftools_output/")


#40kb
wcFST_BtandNonBt_unfilt <- read.table("BtandNonBt_Hzea_variantsonly_wcST_40kb_all.smoothed", header = F)

wcFST_BtandNonBt <- subset(wcFST_BtandNonBt_unfilt, V4 > 10)

#getting significance threshold based on ztransformation
wcFST_BtandNonBt$ztrans <- scale(wcFST_BtandNonBt$V5, center = TRUE, scale = TRUE) #ztransformation 
wcFST_BtandNonBt$SnpName <- seq(1:NROW(wcFST_BtandNonBt))

hist(wcFST_BtandNonBt$ztrans, breaks = 100)

hi_wcFST6 <- subset(wcFST_BtandNonBt, ztrans > 6) #subsetting out any window with an FST greater than 6 SDs from mean.
print(min(hi_wcFST6$V5))
print(mean(wcFST_BtandNonBt$V5))
print(median(wcFST_BtandNonBt$V5))

forPlot_BtandNonBt <- data.frame(cbind(wcFST_BtandNonBt$SnpName,wcFST_BtandNonBt$V1, wcFST_BtandNonBt$V2, wcFST_BtandNonBt$V5))


CMplot(forPlot_BtandNonBt, plot.type="c", r=1.6, cir.legend=TRUE, col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST6$V5)),
       outward=TRUE, cir.legend.col="black", cir.chr.h=.2 ,chr.den.col="pink", file="jpg",
       memo="BtandNonBt_40kb_wcFST", dpi=300, cex.lab = 4, threshold = min(hi_wcFST6$V5), LOG10 = F, chr.labels=seq(1,2958), 
       highlight = hi_wcFST6$SnpName)

#10kb
wcFST_BtandNonBt_unfilt <- read.table("BtandNonBt_Hzea_variantsonly_wcST_10kb_all.smoothed", header = F)

wcFST_BtandNonBt <- subset(wcFST_BtandNonBt_unfilt, V4 > 10)

#getting significance threshold based on ztransformation
wcFST_BtandNonBt$ztrans <- scale(wcFST_BtandNonBt$V5, center = TRUE, scale = TRUE) #ztransformation 
wcFST_BtandNonBt$SnpName <- seq(1:NROW(wcFST_BtandNonBt))

hist(wcFST_BtandNonBt$ztrans, breaks = 100)

hi_wcFST7 <- subset(wcFST_BtandNonBt, ztrans > 6) #subsetting out any window with an FST greater than 6 SDs from mean.
print(min(hi_wcFST7$V5))
print(mean(wcFST_BtandNonBt$V5))
print(median(wcFST_BtandNonBt$V5))

forPlot_BtandNonBt <- data.frame(cbind(wcFST_BtandNonBt$SnpName,wcFST_BtandNonBt$V1, wcFST_BtandNonBt$V2, wcFST_BtandNonBt$V5))


CMplot(forPlot_BtandNonBt, plot.type="c", r=1.6, cir.legend=TRUE, col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST7$V5)),
       outward=TRUE, cir.legend.col="black", cir.chr.h=.2 ,chr.den.col="pink", file="jpg",
       memo="BtandNonBt_10kb_wcFST", dpi=300, cex.lab = 4, threshold = min(hi_wcFST7$V5), LOG10 = F, chr.labels=seq(1,2958), 
       highlight = hi_wcFST7$SnpName)

