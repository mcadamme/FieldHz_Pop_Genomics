#this is my script to analyze my H. zea ddRADseq polymorphism data. 
#M.Fritz 11/23/18

library(LEA);library(OutFLANK);library(adegenet);library(vcfR)

setwd("/media/megan/New Volume/Hz_PopGen_ddRAD_demult/Bowtie_genome_alignments/mpileupANDvcftools_output")

#converting to the proper file format
vcf2geno(input.file = "thinned_FieldHzea2002and2007.recode.vcf", output.file = "FieldHzea2002and2007.geno")
vcf2geno(input.file = "thinned_FieldHzea2007and2012.recode.vcf", output.file = "FieldHzea2007and2012.geno")
vcf2geno(input.file = "thinned_FieldHzea2012and2016.recode.vcf", output.file = "FieldHzea2012and2016.geno")
vcf2geno(input.file = "thinned_FieldHzea2002and2016.recode.vcf", output.file = "FieldHzea2002and2016.geno")


#loading dataset for 2002 and 2007
Field_Hzea_2002and2007.geno_in <- read.fwf("FieldHzea2002and2007.geno", width=rep(1,141))
Field_Hzea_2002and2007.geno <- t(Field_Hzea_2002and2007.geno_in)
Field_Hzea_2002and2007_meta <- read.table("/home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002and2007.txt")

colnames(Field_Hzea_2002and2007_meta) <- c("samp_names")
Field_Hzea_2002and2007_meta$pop <- regmatches(Field_Hzea_2002and2007_meta$samp_names, regexpr("20[[:digit:]]+", Field_Hzea_2002and2007_meta$samp_names))

#OutFLANK analysis for 2002 and 2007
OF_2002and2007_SNPs <- MakeDiploidFSTMat(Field_Hzea_2002and2007.geno, locusNames=seq(1, 27590, by=1), popNames=Field_Hzea_2002and2007_meta$pop) # Provide the transposed genotype file, locus names, and population names.
OF_out1 <- OutFLANK(FstDataFrame=OF_2002and2007_SNPs, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=2, qthreshold=0.05)
OutFLANKResultsPlotter(OF_out1, withOutliers=T, NoCorr=T, Hmin=0.1, binwidth=0.005, Zoom = TRUE, titletext="Scan for selective sweeps 2002 & 2007")

outliers1 <- which(OF_out1$results$OutlierFlag=="TRUE")
print(outliers1)

#getting positions of the outliers
vcf2002and2007 <- read.vcfR("thinned_FieldHzea2002and2007.recode.vcf")
vcfann <- as.data.frame(getFIX(vcf2002and2007))
outliers_2002and2007 <- vcfann[outliers1,]

#write.table(outliers_2002and2007, file = "outliers_2002and2007.txt", col.names = T, row.names = F)

nrow(outliers_2002and2007)
length(unique(outliers_2002and2007$CHROM))
outliers_2002and2007$MARK <- paste(outliers_2002and2007$CHROM, outliers_2002and2007$POS, sep = "_")


#loading dataset for 2007 and 2012
Field_Hzea_2007and2012.geno_in <- read.fwf("FieldHzea2007and2012.geno", width=rep(1,143))
Field_Hzea_2007and2012.geno <- t(Field_Hzea_2007and2012.geno_in)
Field_Hzea_2007and2012_meta <- read.table("/home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007and2012.txt")

colnames(Field_Hzea_2007and2012_meta) <- c("samp_names")
Field_Hzea_2007and2012_meta$pop <- regmatches(Field_Hzea_2007and2012_meta$samp_names, regexpr("20[[:digit:]]+", Field_Hzea_2007and2012_meta$samp_names))

#OutFLANK analysis for 2007 and 2012
OF_2007and2012_SNPs <- MakeDiploidFSTMat(Field_Hzea_2007and2012.geno, locusNames=seq(1, 27590, by=1), popNames=Field_Hzea_2007and2012_meta$pop) 
OF_out2 <- OutFLANK(FstDataFrame=OF_2007and2012_SNPs, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=2, qthreshold=0.05)
OutFLANKResultsPlotter(OF_out2, withOutliers=T, NoCorr=T, Hmin=0.1, binwidth=0.005, Zoom = TRUE, titletext="Scan for selective sweeps 2007 & 2012")

outliers2 <- which(OF_out2$results$OutlierFlag=="TRUE")
print(outliers2)

vcf2007and2012 <- read.vcfR("thinned_FieldHzea2007and2012.recode.vcf")
vcfann <- as.data.frame(getFIX(vcf2007and2012))
outliers_2007and2012 <- vcfann[outliers2,]

#write.table(outliers_2007and2012, file = "outliers_2007and2012.txt", col.names = T, row.names = F)

nrow(outliers_2007and2012)
length(unique(outliers_2007and2012$CHROM))
outliers_2007and2012$MARK <- paste(outliers_2007and2012$CHROM, outliers_2007and2012$POS, sep = "_")



#2012 through 2016
Field_Hzea_2012and2016.geno_in <- read.fwf("FieldHzea2012and2016.geno", width=rep(1,118))
Field_Hzea_2012and2016.geno <- t(Field_Hzea_2012and2016.geno_in)
Field_Hzea_2012and2016_meta <- read.table("/home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012and2016.txt")

colnames(Field_Hzea_2012and2016_meta) <- c("samp_names")
Field_Hzea_2012and2016_meta$pop <- regmatches(Field_Hzea_2012and2016_meta$samp_names, regexpr("20[[:digit:]]+", Field_Hzea_2012and2016_meta$samp_names))

#OutFLANK analysis for 2012 and 2016
OF_2012and2016_SNPs <- MakeDiploidFSTMat(Field_Hzea_2012and2016.geno, locusNames=seq(1, 27590, by=1), popNames=Field_Hzea_2012and2016_meta$pop) 
OF_out3 <- OutFLANK(FstDataFrame=OF_2012and2016_SNPs, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=2, qthreshold=0.05)
OutFLANKResultsPlotter(OF_out3, withOutliers=T, NoCorr=T, Hmin=0.1, binwidth=0.005, Zoom = TRUE, titletext="Scan for selective sweeps 2012 & 2016")

outliers3 <- which(OF_out3$results$OutlierFlag=="TRUE")
print(outliers3)

vcf2012and2016 <- read.vcfR("thinned_FieldHzea2012and2016.recode.vcf")
vcfann <- as.data.frame(getFIX(vcf2012and2016))
outliers_2012and2016 <- vcfann[outliers3,]

#write.table(outliers_2012and2016, file = "outliers_2012and2016.txt", col.names = T, row.names = F)

nrow(outliers_2012and2016)
length(unique(outliers_2012and2016$CHROM))
outliers_2012and2016$MARK <- paste(outliers_2012and2016$CHROM, outliers_2012and2016$POS, sep = "_")



#2002 through 2016
Field_Hzea_2002and2016.geno_in <- read.fwf("FieldHzea2002and2016.geno", width=rep(1,116))
Field_Hzea_2002and2016.geno <- t(Field_Hzea_2002and2016.geno_in)
Field_Hzea_2002and2016_meta <- read.table("/home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002and2016.txt")

colnames(Field_Hzea_2002and2016_meta) <- c("samp_names")
Field_Hzea_2002and2016_meta$pop <- regmatches(Field_Hzea_2002and2016_meta$samp_names, regexpr("20[[:digit:]]+", Field_Hzea_2002and2016_meta$samp_names))

#OutFLANK analysis for 2002 and 2016
OF_2002and2016_SNPs <- MakeDiploidFSTMat(Field_Hzea_2002and2016.geno, locusNames=seq(1, 27590, by=1), popNames=Field_Hzea_2002and2016_meta$pop) 
OF_out4 <- OutFLANK(FstDataFrame=OF_2002and2016_SNPs, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=2, qthreshold=0.05)
OutFLANKResultsPlotter(OF_out4, withOutliers=T, NoCorr=T, Hmin=0.1, binwidth=0.005, Zoom = TRUE, titletext="Scan for selective sweeps 2002 & 2016")

outliers4 <- which(OF_out4$results$OutlierFlag=="TRUE")
print(outliers4)

vcf2002and2016 <- read.vcfR("thinned_FieldHzea2002and2016.recode.vcf")
vcfann <- as.data.frame(getFIX(vcf2002and2016))
outliers_2002and2016 <- vcfann[outliers4,]

#write.table(outliers_2002and2016, file = "outliers_2002and2016.txt", col.names = T, row.names = F)

nrow(outliers_2002and2016)
length(unique(outliers_2002and2016$CHROM))
outliers_2002and2016$MARK <- paste(outliers_2002and2016$CHROM, outliers_2002and2016$POS, sep = "_")

#Outlier Analysis Figures for Pub

#2002&2016 qthreshold 0.05
png("Outflank_Analysis_2002and2016_q0.05.png", units = "px", height = 600, width = 800)
par(mar = c(5,6,4,2))
plot(OF_out4$results$He, OF_out4$results$FST, pch=20, col="grey", ylab = "FST", xlab = "Expected Heterozygosity", 
     ylim = c(0, 0.5), cex.axis = 2, cex.lab = 2.2)
points(OF_out4$results$He[OF_out4$results$qvalues<0.05], OF_out4$results$FST[OF_out4$results$qvalues<0.05], pch=21, col="black")

dev.off()

#qthreshold 0.05
png("Outflank_Analysis_multipanel_years_q0.05.png", units = "px", height = 1200, width = 800)

par(mfrow = c(3,1))
par(mar = c(5,6,4,2))

plot(OF_out1$results$He, OF_out1$results$FST, pch=20, col="grey", ylab = "FST", xlab = "Expected Heterozygosity", 
     ylim = c(0, 0.45), cex.axis = 2, cex.lab = 2.2)
points(OF_out1$results$He[OF_out1$results$qvalues<0.05], OF_out1$results$FST[OF_out1$results$qvalues<0.05], pch=21, col="black")
title("A", cex.main = 3, adj  = 0.05, line = -5)

plot(OF_out2$results$He, OF_out2$results$FST, pch=20, col="grey", ylab = "FST", xlab = "Expected Heterozygosity", 
     ylim = c(0, 0.45), cex.axis = 2, cex.lab = 2.2)
points(OF_out2$results$He[OF_out2$results$qvalues<0.05], OF_out2$results$FST[OF_out2$results$qvalues<0.05], pch=21, col="black")
title("B", cex.main = 3, adj  = 0.05, line = -5)

plot(OF_out3$results$He, OF_out3$results$FST, pch=20, col="grey", ylab = "FST", xlab = "Expected Heterozygosity", 
     ylim = c(0, 0.45), cex.axis = 2, cex.lab = 2.2)
points(OF_out3$results$He[OF_out3$results$qvalues<0.05], OF_out3$results$FST[OF_out3$results$qvalues<0.05], pch=21, col="black")
title("C", cex.main = 3, adj  = 0.05, line = -5)

dev.off()

#qthreshold 0.01
png("Outflank_Analysis_multipanel_years_q0.01.png", units = "px", height = 1200, width = 800)

par(mfrow = c(3,1))
par(mar = c(5,6,4,2))

plot(OF_out1$results$He, OF_out1$results$FST, pch=20, col="grey", ylab = "FST", xlab = "Expected Heterozygosity", 
     ylim = c(0, 0.45), cex.axis = 2, cex.lab = 2.2)
points(OF_out1$results$He[OF_out1$results$qvalues<0.01], OF_out1$results$FST[OF_out1$results$qvalues<0.01], pch=21, col="black")
title("A", cex.main = 3, adj  = 0.05, line = -5)

plot(OF_out2$results$He, OF_out2$results$FST, pch=20, col="grey", ylab = "FST", xlab = "Expected Heterozygosity", 
     ylim = c(0, 0.45), cex.axis = 2, cex.lab = 2.2)
points(OF_out2$results$He[OF_out2$results$qvalues<0.01], OF_out2$results$FST[OF_out2$results$qvalues<0.01], pch=21, col="black")
title("B", cex.main = 3, adj  = 0.05, line = -5)

plot(OF_out3$results$He, OF_out3$results$FST, pch=20, col="grey", ylab = "FST", xlab = "Expected Heterozygosity", 
     ylim = c(0, 0.45), cex.axis = 2, cex.lab = 2.2)
points(OF_out3$results$He[OF_out3$results$qvalues<0.01], OF_out3$results$FST[OF_out3$results$qvalues<0.01], pch=21, col="black")
title("C", cex.main = 3, adj  = 0.05, line = -5)

dev.off()


#subsetting by top sig outliers in two by-year comparisons

outliers_firsttwo <- merge(outliers_2002and2007,outliers_2002and2016, by = "MARK")
length(unique(outliers_firsttwo$MARK))
print(unique(outliers_firsttwo$MARK))

outliers_nexttwo <- merge(outliers_2007and2012,outliers_2002and2016, by = "MARK")
length(unique(outliers_nexttwo$MARK))
print(unique(outliers_nexttwo$MARK))

outliers_lasttwo <- merge(outliers_2012and2016,outliers_2002and2016, by = "MARK")
length(unique(outliers_lasttwo$MARK))
print(unique(outliers_lasttwo$MARK))

##looking at shared outliers across mult by-year comparisons
outliers_all1 <- merge(outliers_firsttwo,outliers_nexttwo, by = "MARK")
length(unique(outliers_all1$MARK))
print(unique(outliers_all1$MARK))

outliers_all2 <- merge(outliers_lasttwo,outliers_nexttwo, by = "MARK")
length(unique(outliers_all2$MARK))
print(unique(outliers_all2$MARK))
#this outlier has weird changes in allele frequency over time, inconsistent with a selective sweep.


