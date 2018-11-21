#this is my script to analyze my H. zea ddRADseq polymorphism data. 


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
nrow(outliers_2002and2007)
length(unique(outliers_2002and2007$CHROM))

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
nrow(outliers_2007and2012)
length(unique(outliers_2007and2012$CHROM))


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
nrow(outliers_2012and2016)
length(unique(outliers_2012and2016$CHROM))


#2012 through 2016
Field_Hzea_2012and2016.geno_in <- read.fwf("FieldHzea2012and2016.geno", width=rep(1,118))
Field_Hzea_2012and2016.geno <- t(Field_Hzea_2012and2016.geno_in)
Field_Hzea_2012and2016_meta <- read.table("/home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012and2016.txt")

colnames(Field_Hzea_2012and2016_meta) <- c("samp_names")
Field_Hzea_2012and2016_meta$pop <- regmatches(Field_Hzea_2012and2016_meta$samp_names, regexpr("20[[:digit:]]+", Field_Hzea_2012and2016_meta$samp_names))

#OutFLANK analysis for 2002 and 2016
OF_2002and2016_SNPs <- MakeDiploidFSTMat(Field_Hzea_2002and2016.geno, locusNames=seq(1, 27590, by=1), popNames=Field_Hzea_2002and2016_meta$pop) 
OF_out4 <- OutFLANK(FstDataFrame=OF_2002and2016_SNPs, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=2, qthreshold=0.05)
OutFLANKResultsPlotter(OF_out4, withOutliers=T, NoCorr=T, Hmin=0.1, binwidth=0.005, Zoom = TRUE, titletext="Scan for selective sweeps 2012 & 2016")

outliers3 <- which(OF_out3$results$OutlierFlag=="TRUE")
print(outliers3)

vcf2012and2016 <- read.vcfR("thinned_FieldHzea2012and2016.recode.vcf")
vcfann <- as.data.frame(getFIX(vcf2012and2016))
outliers_2012and2016 <- vcfann[outliers3,]
nrow(outliers_2012and2016)
length(unique(outliers_2012and2016$CHROM))


#Outlier Analysis Figure for Pub

#qthreshold 0.05
png("Outflank_Analysis_multipanel_years_q0.05.png", units = "px", height = 1200, width = 800)

par(mfrow = c(3,1))

plot(OF_out1$results$He, OF_out1$results$FST, pch=20, col="grey", ylab = "FST", xlab = "Expected Heterozygosity", 
     ylim = c(0, 0.45), cex.axis = 1.5, cex.lab = 1.2)
points(OF_out1$results$He[OF_out1$results$qvalues<0.05], OF_out1$results$FST[OF_out1$results$qvalues<0.05], pch=21, col="black")
title("A", cex.main = 1.5, adj  = 0.05, line = -3)

plot(OF_out2$results$He, OF_out2$results$FST, pch=20, col="grey", ylab = "FST", xlab = "Expected Heterozygosity", 
     ylim = c(0, 0.45), cex.axis = 1.5, cex.lab = 1.2)
points(OF_out2$results$He[OF_out2$results$qvalues<0.05], OF_out2$results$FST[OF_out2$results$qvalues<0.05], pch=21, col="black")
title("B", cex.main = 1.5, adj  = 0.05, line = -3)

plot(OF_out3$results$He, OF_out3$results$FST, pch=20, col="grey", ylab = "FST", xlab = "Expected Heterozygosity", 
     ylim = c(0, 0.45), cex.axis = 1.5, cex.lab = 1.2)
points(OF_out3$results$He[OF_out3$results$qvalues<0.05], OF_out3$results$FST[OF_out3$results$qvalues<0.05], pch=21, col="black")
title("C", cex.main = 1.5, adj  = 0.05, line = -3)

dev.off()

#qthreshold 0.01
png("Outflank_Analysis_multipanel_years_q0.01.png", units = "px", height = 1200, width = 800)

par(mfrow = c(3,1))
plot(OF_out1$results$He, OF_out1$results$FST, pch=20, col="grey", ylab = "FST", xlab = "Expected Heterozygosity", 
     ylim = c(0, 0.45), cex.axis = 1.5, cex.lab = 1.2)
points(OF_out1$results$He[OF_out1$results$qvalues<0.01], OF_out1$results$FST[OF_out1$results$qvalues<0.01], pch=21, col="black")
title("A", cex.main = 1.5, adj  = 0.05, line = -3)

plot(OF_out2$results$He, OF_out2$results$FST, pch=20, col="grey", ylab = "FST", xlab = "Expected Heterozygosity", 
     ylim = c(0, 0.45), cex.axis = 1.5, cex.lab = 1.2)
points(OF_out2$results$He[OF_out2$results$qvalues<0.01], OF_out2$results$FST[OF_out2$results$qvalues<0.01], pch=21, col="black")
title("B", cex.main = 1.5, adj  = 0.05, line = -3)

plot(OF_out3$results$He, OF_out3$results$FST, pch=20, col="grey", ylab = "FST", xlab = "Expected Heterozygosity", 
     ylim = c(0, 0.45), cex.axis = 1.5, cex.lab = 1.2)
points(OF_out3$results$He[OF_out3$results$qvalues<0.01], OF_out3$results$FST[OF_out3$results$qvalues<0.01], pch=21, col="black")
title("C", cex.main = 1.5, adj  = 0.05, line = -3)

dev.off()


#subsetting by top sig outliers
hioutliers3 <- which(OF_out3$results$OutlierFlag=="TRUE" & OF_out3$results$FST > 0.15)
print(hioutliers3)
hi_outliers_2002and2016 <- vcfann[hioutliers3,]
length(unique(hi_outliers_2002and2016$CHROM))

outliers_firsttwo <- merge(outliers_2002and2007,outliers_2002and2012, by = "CHROM")
length(unique(outliers_firsttwo$CHROM))
print(unique(outliers_firsttwo$CHROM))

outliers_lasttwo <- merge(outliers_2002and2012,outliers_2002and2016, by = "CHROM")
length(unique(outliers_lasttwo$CHROM))
print(unique(outliers_lasttwo$CHROM))

outliers_all <- merge(outliers_firsttwo,outliers_lasttwo, by = "CHROM")
length(unique(outliers_all$CHROM))
print(unique(outliers_all$CHROM))


