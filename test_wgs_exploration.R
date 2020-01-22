
library(LEA);library(OutFLANK);library(adegenet);library(vcfR);library(CMplot)


setwd("/media/megan/New Volume/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1")


#first 2002-2017 comparison
#40kb
wcFST_2002_2017 <- read.table("2002and2017_40kb_wcFST_all.smoothed", header = F)

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
       memo="", dpi=300, cex.lab = 4, threshold = min(hi_wcFST0$V5), LOG10 = F, chr.labels=seq(1,2958), 
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
       memo="", dpi=300, cex.lab = 4, threshold = min(hi_wcFST1$V5), LOG10 = F, chr.labels=seq(1,2958), 
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
       memo="", dpi=300, cex.lab = 4, threshold = min(hi_wcFST2$V5), LOG10 = F, chr.labels=seq(1,2958), 
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
       memo="", dpi=300, cex.lab = 4, threshold = min(hi_wcFST3$V5), LOG10 = F, chr.labels=seq(1,2958), 
       highlight = hi_wcFST3$SnpName)



#2002-2012 comparison
#40kb
wcFST_2002_2012 <- read.table("2002and2012_40kb_wcFST_all.smoothed", header = F)

#getting significance threshold based on ztransformation
wcFST_2002_2012$ztrans <- scale(wcFST_2002_2012$V5, center = TRUE, scale = TRUE) #ztransformation 
wcFST_2002_2012$SnpName <- seq(1:NROW(wcFST_2002_2012))

hist(wcFST_2002_2012$ztrans, breaks = 100)

hi_wcFST4 <- subset(wcFST_2002_2012, ztrans > 6) #subsetting out any window with an FST greater than 6 SDs from mean.
print(min(hi_wcFST4$V5))


forPlot_2002_2012 <- data.frame(cbind(wcFST_2002_2012$SnpName,wcFST_2002_2012$V1, wcFST_2002_2012$V2, wcFST_2002_2012$V5))

#2002-2012 Circular Manhattan plot function
CMplot(forPlot_2002_2012, plot.type="c", r=1.6, cir.legend=TRUE, col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST4$V5)),
       outward=TRUE, cir.legend.col="black", cir.chr.h=.2 ,chr.den.col="pink", file="jpg",
       memo="", dpi=300, cex.lab = 4, threshold = min(hi_wcFST4$V5), LOG10 = F, chr.labels=seq(1,2958), highlight = hi_wcFST4$SnpName)


#2012-2017 comparison
#40kb
wcFST_2012_2017 <- read.table("2012and2017_40kb_wcFST_all.smoothed", header = F)

#getting significance threshold based on ztransformation
wcFST_2012_2017$ztrans <- scale(wcFST_2012_2017$V5, center = TRUE, scale = TRUE) #ztransformation 
wcFST_2012_2017$SnpName <- seq(1:NROW(wcFST_2012_2017))

hist(wcFST_2012_2017$ztrans, breaks = 100)

hi_wcFST5 <- subset(wcFST_2012_2017, ztrans > 6) #subsetting out any window with an FST greater than 6 SDs from mean.
print(min(hi_wcFST5$V5))


forPlot_2012_2017 <- data.frame(cbind(wcFST_2012_2017$SnpName,wcFST_2012_2017$V1, wcFST_2012_2017$V2, wcFST_2012_2017$V5))

#2002-2012 Circular Manhattan plot function
CMplot(forPlot_2012_2017, plot.type="c", r=1.6, cir.legend=TRUE, col = c("grey30", "grey60"), cex = 0.8, ylim = c(0,max(hi_wcFST5$V5)),
       outward=TRUE, cir.legend.col="black", cir.chr.h=.2 ,chr.den.col="pink", file="jpg",
       memo="", dpi=300, cex.lab = 4, threshold = min(hi_wcFST5$V5), LOG10 = F, chr.labels=seq(1,2958), highlight = hi_wcFST5$SnpName)



#Loading scafname and gff files for getting genes nearby outliers

#loading identifier file
scaf_names <- read.table("/media/megan/New Volume/Hz_PopGen_ddRAD_demult/Bowtie_genome_alignments/Hzea_genome/scaffold_and_contig_names.txt",
                         header = T)



#Now pulling out the parts of the gff3 file that are important
gff3 <- read.table("/media/megan/New Volume/Hz_PopGen_ddRAD_demult/Bowtie_genome_alignments/Hzea_genome/HzOGS2-15205-fixed_note-added.gff3",
                   sep="\t", stringsAsFactors = F)
head(gff3)

outliers_allpops_withIDs <- merge(hi_wcFST, scaf_names, by.x = "V1", by.y = "CHROM")
outliers_allyears_withAnnot <- merge(outliers_allpops_withIDs, gff3, by.x = "Scaf_ID", by.y = "V1")

print(unique(outliers_allyears_withAnnot$V9))
cad86c <- subset(outliers_allyears_withAnnot, V3.y == "gene" & Scaf_ID == "scaffold_20" & V2.x > 500000 & V2.x < 700000)
sub_outliers_allyears_withAnnot <- subset(outliers_allyears_withAnnot, V3.y == "gene" & V5.y > (V2.x -50000) & V4.y < (V2.x + 50000))

KZ118005 <- subset(outliers_allyears_withAnnot, V3.y == "gene" & V1 == "KZ118005.1" & V2.x >165000 & V2.x < 190000)
print(unique(KZ118005$V9))
