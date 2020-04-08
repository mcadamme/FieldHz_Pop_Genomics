#this is my script to analyze my H. zea ddRADseq polymorphism data. 
#M.Fritz 11/23/18

library(LEA);library(OutFLANK);library(adegenet);library(vcfR)



#Loading scafname and gff files for getting genes nearby outliers

#loading identifier file
scaf_names <- read.table("/media/megan/New Volume1/Hz_PopGen_ddRAD_demult/Bowtie_genome_alignments/Hzea_genome/scaffold_and_contig_names.txt",
                         header = T)



#Now pulling out the parts of the gff3 file that are important
gff3 <- read.table("/media/megan/New Volume1/Hz_PopGen_ddRAD_demult/Bowtie_genome_alignments/Hzea_genome/HzOGS2-15205-fixed_note-added.gff3",
                   sep="\t", stringsAsFactors = F)
head(gff3)


#Outlier Analysis of BCFtools output

setwd("/media/megan/New Volume1/Hz_PopGen_ddRAD_demult/Bowtie_genome_alignments/mpileupANDvcftools_output")

#converting to the proper file format
#overall pop across years
vcf2geno(input.file = "thinned_FieldHzea_allpops.recode.vcf", output.file = "FieldHzea_allpops.geno")

Field_Hzea_allpops.geno_in <- read.fwf("FieldHzea_allpops.geno", width=rep(1,259))
Field_Hzea_allpops.geno <- t(Field_Hzea_allpops.geno_in)
Field_Hzea_allpops_meta <- read.table("/home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/allpops_no_low.txt")

colnames(Field_Hzea_allpops_meta) <- c("samp_names")
Field_Hzea_allpops_meta$pop <- regmatches(Field_Hzea_allpops_meta$samp_names, regexpr("20[[:digit:]]+", Field_Hzea_allpops_meta$samp_names))

#OutFLANK analysis for all years
OF_allpops_SNPs <- MakeDiploidFSTMat(Field_Hzea_allpops.geno, locusNames=seq(1, 14398, by=1), popNames=Field_Hzea_allpops_meta$pop) # Provide the transposed genotype file, locus names, and population names.
OF_out <- OutFLANK(FstDataFrame=OF_allpops_SNPs, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=4, qthreshold=0.05)
OutFLANKResultsPlotter(OF_out, withOutliers=T, NoCorr=T, Hmin=0.1, binwidth=0.005, Zoom = TRUE, titletext="Scan for selective sweeps across years")


#getting positions of the outliers
vcf_allpops <- read.vcfR("thinned_FieldHzea_allpops.recode.vcf")
vcfann <- as.data.frame(getFIX(vcf_allpops))

outliers_allyears_df <- data.frame(OF_out$results, header = T)
outliers_allyears_df_withpos <- cbind(vcfann[,c(1:2)],outliers_allyears_df)
outliers_allyears <- subset(outliers_allyears_df_withpos, OutlierFlag == "TRUE")

outliers_allyears_uniqScafs <- unique(outliers_allyears$CHROM)
outliers_allyears_uniqScafs

write.table(outliers_allyears, file = "FileS1_outliers_allyears.txt", col.names = T, row.names = F)

min(outliers_allyears$FST)
max(outliers_allyears$FST)

#all_years qval 0.05
png("Fig2_Outflank_Analysis_allyears_qval0.05.png", units = "px", height = 600, width = 800)
par(mar = c(5,6,4,2))
plot(OF_out$results$He, OF_out$results$FST, pch=20, col="grey", ylab = "FST", xlab = "Expected Heterozygosity", 
     ylim = c(0, 0.3), cex.axis = 2, cex.lab = 2.2)
points(OF_out$results$He[OF_out$results$qvalues<0.05], OF_out$results$FST[OF_out$results$qvalues<0.05], pch=21, col="black")

dev.off()

#Analysis by pair of years - just for curiosity.  Did not go into paper.
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
OF_2002and2007_SNPs <- MakeDiploidFSTMat(Field_Hzea_2002and2007.geno, locusNames=seq(1, 14398, by=1), popNames=Field_Hzea_2002and2007_meta$pop) # Provide the transposed genotype file, locus names, and population names.
OF_out1 <- OutFLANK(FstDataFrame=OF_2002and2007_SNPs, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=2, qthreshold=0.1)
OutFLANKResultsPlotter(OF_out1, withOutliers=T, NoCorr=T, Hmin=0.1, binwidth=0.005, Zoom = TRUE, titletext="Scan for selective sweeps 2002 & 2007")

outliers1 <- which(OF_out1$results$OutlierFlag=="TRUE")
print(outliers1)

#getting positions of the outliers
vcf2002and2007 <- read.vcfR("thinned_FieldHzea2002and2007.recode.vcf")
vcfann <- as.data.frame(getFIX(vcf2002and2007))
outliers_2002and2007 <- vcfann[outliers1,]

write.table(outliers_2002and2007, file = "outliers_2002and2007.txt", col.names = T, row.names = F)

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
OF_2007and2012_SNPs <- MakeDiploidFSTMat(Field_Hzea_2007and2012.geno, locusNames=seq(1, 14398, by=1), popNames=Field_Hzea_2007and2012_meta$pop) 
OF_out2 <- OutFLANK(FstDataFrame=OF_2007and2012_SNPs, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=2, qthreshold=0.1)
OutFLANKResultsPlotter(OF_out2, withOutliers=T, NoCorr=T, Hmin=0.1, binwidth=0.005, Zoom = TRUE, titletext="Scan for selective sweeps 2007 & 2012")

outliers2 <- which(OF_out2$results$OutlierFlag=="TRUE")
print(outliers2)

vcf2007and2012 <- read.vcfR("thinned_FieldHzea2007and2012.recode.vcf")
vcfann <- as.data.frame(getFIX(vcf2007and2012))
outliers_2007and2012 <- vcfann[outliers2,]

write.table(outliers_2007and2012, file = "outliers_2007and2012.txt", col.names = T, row.names = F)

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
OF_2012and2016_SNPs <- MakeDiploidFSTMat(Field_Hzea_2012and2016.geno, locusNames=seq(1, 14398, by=1), popNames=Field_Hzea_2012and2016_meta$pop) 
OF_out3 <- OutFLANK(FstDataFrame=OF_2012and2016_SNPs, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=2, qthreshold=0.1)
OutFLANKResultsPlotter(OF_out3, withOutliers=T, NoCorr=T, Hmin=0.1, binwidth=0.005, Zoom = TRUE, titletext="Scan for selective sweeps 2012 & 2016")

outliers3 <- which(OF_out3$results$OutlierFlag=="TRUE")
print(outliers3)

vcf2012and2016 <- read.vcfR("thinned_FieldHzea2012and2016.recode.vcf")
vcfann <- as.data.frame(getFIX(vcf2012and2016))
outliers_2012and2016 <- vcfann[outliers3,]

write.table(outliers_2012and2016, file = "outliers_2012and2016.txt", col.names = T, row.names = F)

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
OF_2002and2016_SNPs <- MakeDiploidFSTMat(Field_Hzea_2002and2016.geno, locusNames=seq(1, 14398, by=1), popNames=Field_Hzea_2002and2016_meta$pop) 
OF_out4 <- OutFLANK(FstDataFrame=OF_2002and2016_SNPs, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=2, qthreshold=0.1)
OutFLANKResultsPlotter(OF_out4, withOutliers=T, NoCorr=T, Hmin=0.1, binwidth=0.005, Zoom = TRUE, titletext="Scan for selective sweeps 2002 & 2016")

outliers4 <- which(OF_out4$results$OutlierFlag=="TRUE")
print(outliers4)
outliers4_alpha0.01 <-which(OF_out4$results$q < 0.01)

vcf2002and2016 <- read.vcfR("thinned_FieldHzea2002and2016.recode.vcf")
vcfann <- as.data.frame(getFIX(vcf2002and2016))
outliers_2002and2016 <- vcfann[outliers4,]
outliers_2002and2016_q0.01 <- vcfann[outliers4_alpha0.01,]

write.table(outliers_2002and2016, file = "outliers_2002and2016.txt", col.names = T, row.names = F)

nrow(outliers_2002and2016)
length(unique(outliers_2002and2016$CHROM))
outliers_2002and2016$MARK <- paste(outliers_2002and2016$CHROM, outliers_2002and2016$POS, sep = "_")

#Outlier Analysis Figures for Pub

#2002&2016 qval 0.05
png("Fig2_Outflank_Analysis_2002and2016_qval0.05.png", units = "px", height = 600, width = 800)
par(mar = c(5,6,4,2))
plot(OF_out4$results$He, OF_out4$results$FST, pch=20, col="grey", ylab = "FST", xlab = "Expected Heterozygosity", 
     ylim = c(0, 0.5), cex.axis = 2, cex.lab = 2.2)
points(OF_out4$results$He[OF_out4$results$qvalues<0.05], OF_out4$results$FST[OF_out4$results$qvalues<0.05], pch=21, col="black")

dev.off()

#qthreshold 0.05
png("Outflank_Analysis_multipanel_years_qval0.05.png", units = "px", height = 1200, width = 800)

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
png("Outflank_Analysis_multipanel_years_qval0.01.png", units = "px", height = 1200, width = 800)

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

#adding scaffold names to outlier dataframes for getting nearby genes
outliers_allpops_withIDs <- merge(outliers_allpops, scaf_names, by = "CHROM")
outliers_2002and2007_withIDs <- merge(outliers_2002and2007, scaf_names, by = "CHROM")
outliers_2007and2012_withIDs <- merge(outliers_2007and2012, scaf_names, by = "CHROM")
outliers_2012and2016_withIDs <- merge(outliers_2012and2016, scaf_names, by = "CHROM")
outliers_2002and2016_withIDs <- merge(outliers_2002and2016, scaf_names, by = "CHROM")

#all years
outliers_allyears_withAnnot <- merge(outliers_allpops_withIDs, gff3, by.x = "Scaf_ID", by.y = "V1")
outliers_allyears_withAnnot$POS <- as.integer(as.character(outliers_allyears_withAnnot$POS))
outliers_allyears_withAnnot <- subset(outliers_allyears_withAnnot, V3 == "gene" & V5 > (POS -10000) & V4 < (POS + 10000))
write.table(outliers_allyears_withAnnot, file = "Sig_outliers_withAnnot_allyears.txt", row.names = F)
write.table(unique(outliers_allyears_withAnnot$V9), file = "Sig_outliers_withAnnot_allyears_genesOnly.txt", row.names = F)


#2002 and 2007
outliers_2002and2007_withAnnot <- merge(outliers_2002and2007_withIDs, gff3, by.x = "Scaf_ID", by.y = "V1")
outliers_2002and2007_withAnnot$POS <- as.integer(as.character(outliers_2002and2007_withAnnot$POS))
outliers_2002and2007_withAnnot <- subset(outliers_2002and2007_withAnnot, V3 == "gene" & V5 > (POS -10000) & V4 < (POS + 10000))
write.table(outliers_2002and2007_withAnnot, file = "Sig_outliers_withAnnot_2002and2007.txt", row.names = F)
write.table(unique(outliers_2002and2007_withAnnot$V9), file = "Sig_outliers_withAnnot_2002and2007_genesOnly.txt", row.names = F)

#2007 and 2012
outliers_2007and2012_withAnnot <- merge(outliers_2007and2012_withIDs, gff3, by.x = "Scaf_ID", by.y = "V1")
outliers_2007and2012_withAnnot$POS <- as.integer(as.character(outliers_2007and2012_withAnnot$POS))
outliers_2007and2012_withAnnot <- subset(outliers_2007and2012_withAnnot, V3 == "gene" & V5 > (POS -10000) & V4 < (POS + 10000))
write.table(outliers_2007and2012_withAnnot, file = "Sig_outliers_withAnnot_2007and2012.txt", row.names = F)
write.table(unique(outliers_2007and2012_withAnnot$V9), file = "Sig_outliers_withAnnot_2007and2012_genesOnly.txt", row.names = F)


#2012 and 2016
outliers_2012and2016_withAnnot <- merge(outliers_2012and2016_withIDs, gff3, by.x = "Scaf_ID", by.y = "V1")
outliers_2012and2016_withAnnot$POS <- as.integer(as.character(outliers_2012and2016_withAnnot$POS))
outliers_2012and2016_withAnnot <- subset(outliers_2012and2016_withAnnot, V3 == "gene" & V5 > (POS - 10000) & V4 < (POS + 10000))
write.table(outliers_2012and2016_withAnnot, file = "Sig_outliers_withAnnot_2012and2016.txt", row.names = F)
write.table(unique(outliers_2012and2016_withAnnot$V9), file = "Sig_outliers_withAnnot_2012and2016_genesOnly.txt", row.names = F)


#2002 and 2016
outliers_2002and2016_withAnnot <- merge(outliers_2002and2016_withIDs, gff3, by.x = "Scaf_ID", by.y = "V1")
outliers_2002and2016_withAnnot$POS <- as.integer(as.character(outliers_2002and2016_withAnnot$POS))
outliers_2002and2016_withAnnot <- subset(outliers_2002and2016_withAnnot, V3 == "gene" & V5 > (POS -10000) & V4 < (POS + 10000))
write.table(outliers_2002and2016_withAnnot, file = "Sig_outliers_withAnnot_2002and2016.txt", row.names = F)
write.table(unique(outliers_2002and2016_withAnnot$V9), file = "Sig_outliers_withAnnot_2002and2016_genesOnly.txt", row.names = F)
